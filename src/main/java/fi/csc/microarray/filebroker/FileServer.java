package fi.csc.microarray.filebroker;

import java.io.File;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Arrays;
import java.util.Timer;
import java.util.concurrent.TimeUnit;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.manager.ManagerClient;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingListener;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.message.BooleanMessage;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;
import fi.csc.microarray.messaging.message.UrlMessage;
import fi.csc.microarray.service.KeepAliveShutdownHandler;
import fi.csc.microarray.service.ShutdownCallback;
import fi.csc.microarray.util.FileCleanUpTimerTask;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.MemUtil;

public class FileServer extends NodeBase implements MessagingListener, ShutdownCallback {
	/**
	 * Logger for this class
	 */
	private static Logger logger;

	private MessagingEndpoint endpoint;	
	private ManagerClient managerClient;
	private AuthorisedUrlRepository urlRepository;

	private File userDataRoot;
	private String publicDataPath;
	private String host;
	private int port;

	private int cleanUpTriggerLimitPercentage;
	private int cleanUpTargetPercentage;
	private int cleanUpMinimumFileAge;
	private long minimumSpaceForAcceptUpload;

	public static void main(String[] args) {
		// we should be able to specify alternative user dir for testing... and replace maybe that previous hack
		//DirectoryLayout.getInstance(new File("chipster-userdir-fileserver")).getConfiguration();
		
		DirectoryLayout.getInstance().getConfiguration();
		new FileServer(null);
	}
	
    public FileServer(String configURL) {

    	try {
    		// initialise dir and logging
    		DirectoryLayout.initialiseServerLayout(
    		        Arrays.asList(new String[] {"frontend", "filebroker"}),
    		        configURL);
    		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
    		logger = Logger.getLogger(FileServer.class);

    		// initialise url repository
    		File fileRepository = DirectoryLayout.getInstance().getFileRoot();
    		this.host = configuration.getString("filebroker", "url");
    		this.port = configuration.getInt("filebroker", "port");
    		
    		this.urlRepository = new AuthorisedUrlRepository(host, port);
    		this.publicDataPath = configuration.getString("filebroker", "public-data-path");

    		// boot up file server
    		JettyFileServer fileServer = new JettyFileServer(urlRepository);
    		fileServer.start(fileRepository.getPath(), port);

    		// start scheduler
    		String userDataPath = configuration.getString("filebroker", "user-data-path");
    		userDataRoot = new File(fileRepository, userDataPath);
    		cleanUpTriggerLimitPercentage = configuration.getInt("filebroker", "clean-up-trigger-limit-percentage");
    		cleanUpTargetPercentage = configuration.getInt("filebroker", "clean-up-target-percentage");
    		cleanUpMinimumFileAge = configuration.getInt("filebroker", "clean-up-minimum-file-age");
    		minimumSpaceForAcceptUpload = 1024*1024*configuration.getInt("filebroker", "minimum-space-for-accept-upload");
    		
    		int cutoff = 1000 * configuration.getInt("filebroker", "file-life-time");
    		int cleanUpFrequency = 1000 * configuration.getInt("filebroker", "clean-up-frequency");
    		int checkFrequency = 1000 * 5;
    		Timer t = new Timer("frontend-scheduled-tasks", true);
    		t.schedule(new FileCleanUpTimerTask(userDataRoot, cutoff), 0, cleanUpFrequency);
    		t.schedule(new JettyCheckTimerTask(fileServer), 0, checkFrequency);

    		// initialise messaging
    		this.endpoint = new MessagingEndpoint(this);
    		MessagingTopic filebrokerTopic = endpoint.createTopic(Topics.Name.AUTHORISED_FILEBROKER_TOPIC, AccessMode.READ);
    		filebrokerTopic.setListener(this);
    		this.managerClient = new ManagerClient(endpoint); 

    		// create keep-alive thread and register shutdown hook
    		KeepAliveShutdownHandler.init(this);

    		
    		logger.info("fileserver is up and running [" + ApplicationConstants.VERSION + "]");
    		logger.info("[mem: " + MemUtil.getMemInfo() + "]");

    	} catch (Exception e) {
    		e.printStackTrace();
    		logger.error(e, e);
    	}
    }


	public String getName() {
		return "filebroker";
	}


	public void onChipsterMessage(ChipsterMessage msg) {
		try {

			if (msg instanceof CommandMessage && CommandMessage.COMMAND_URL_REQUEST.equals(((CommandMessage)msg).getCommand())) {
				CommandMessage requestMessage = (CommandMessage) msg;
				boolean useCompression = requestMessage.getParameters().contains(ParameterMessage.PARAMETER_USE_COMPRESSION);
				URL url = urlRepository.createAuthorisedUrl(useCompression);
				UrlMessage reply = new UrlMessage(url);
				endpoint.replyToMessage(msg, reply);
				managerClient.urlRequest(msg.getUsername(), url);
				
			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_PUBLIC_URL_REQUEST.equals(((CommandMessage)msg).getCommand())) {
				URL url = getPublicUrL();
				UrlMessage reply = new UrlMessage(url);
				endpoint.replyToMessage(msg, reply);
				managerClient.publicUrlRequest(msg.getUsername(), url);

			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_DISK_SPACE_REQUEST.equals(((CommandMessage)msg).getCommand())) {
				handleSpaceRequest((CommandMessage)msg);

			} else {
				logger.error("message " + msg.getMessageID() + " not understood");
			}
			
		} catch (Exception e) {
			logger.error(e, e);
		}
	}

	private void handleSpaceRequest(CommandMessage requestMessage) throws JMSException {
		long size = Long.parseLong(requestMessage.getNamedParameter(ParameterMessage.PARAMETER_DISK_SPACE));
		logger.debug("disk space request for " + size + " bytes");

		long usableSpaceSoftLimit =  (long) ((double)userDataRoot.getTotalSpace()*(double)(100-cleanUpTriggerLimitPercentage)/100);
		long usableSpaceHardLimit = minimumSpaceForAcceptUpload;
		
		// deal with the weird config case of soft limit being smaller than hard limit
		if (usableSpaceSoftLimit < usableSpaceHardLimit) {
			usableSpaceSoftLimit = usableSpaceHardLimit;
		}
		
		boolean spaceAvailable;
		
		// space available, clean up limit will not be reached
		if (userDataRoot.getUsableSpace() - size >= usableSpaceSoftLimit) {
			logger.debug("enough space available, no need to do anything");
			spaceAvailable = true;
		} 

		// space available, clean up soft limit will be reached, hard will not be reached
		else if (userDataRoot.getUsableSpace() - size >= usableSpaceHardLimit) {
			logger.debug("space available, more preferred");
			spaceAvailable = true;
		
			// schedule clean up
			new Thread(new Runnable() {
				@Override
				public void run() {
					try {
						long cleanUpBeginTime = System.currentTimeMillis();
						Files.makeSpaceInDirectory(userDataRoot, 100-cleanUpTargetPercentage, cleanUpMinimumFileAge, TimeUnit.SECONDS);
						logger.info("cache cleanup took " + (System.currentTimeMillis() - cleanUpBeginTime) + " ms");
					} catch (Exception e) {
						logger.warn("exception while cleaning cache", e);
					}
				}
			}, "chipster-fileserver-cache-cleanup").start();
		} 
		
		// hard limit will be reached, try to make more immediately
		else if (userDataRoot.getUsableSpace() - size > 0){
			logger.debug("making space");
			try {
				long cleanUpBeginTime = System.currentTimeMillis();
				Files.makeSpaceInDirectory(userDataRoot, 100-cleanUpTargetPercentage, cleanUpMinimumFileAge, TimeUnit.SECONDS);
				logger.info("cache cleanup took " + (System.currentTimeMillis() - cleanUpBeginTime) + " ms");
			} catch (Exception e) {
				logger.warn("exception while cleaning cache", e);
			}
			logger.debug("usable after cleaning: " + userDataRoot.getUsableSpace());
			logger.debug("minimum extra: " + minimumSpaceForAcceptUpload);

			// check if cleaned up enough 
			if (userDataRoot.getUsableSpace() >= size + minimumSpaceForAcceptUpload ) {
				logger.debug("enough after cleaning");
				spaceAvailable = true;
			} else {
				logger.debug("not enough after cleaning");
				spaceAvailable = false;
			}
		} 
		
		// request more than total, no can do
		else {
			spaceAvailable = false;
		}
		
		// send reply
		BooleanMessage reply = new BooleanMessage(spaceAvailable);
		endpoint.replyToMessage(requestMessage, reply);
	}

	public void shutdown() {
		logger.info("shutdown requested");

		// close messaging endpoint
		try {
			this.endpoint.close();
		} catch (JMSException e) {
			logger.error("closing messaging endpoint failed", e);
		}

		logger.info("shutting down");
	}

	public URL getPublicUrL() throws MalformedURLException {
		return new URL(host + ":" + port + "/" + publicDataPath);		
	}
}
