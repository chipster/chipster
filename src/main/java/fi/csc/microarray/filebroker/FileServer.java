package fi.csc.microarray.filebroker;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.TimeUnit;

import javax.jms.JMSException;

import org.apache.commons.io.FileUtils;
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
import fi.csc.microarray.messaging.message.UrlListMessage;
import fi.csc.microarray.messaging.message.UrlMessage;
import fi.csc.microarray.service.KeepAliveShutdownHandler;
import fi.csc.microarray.service.ShutdownCallback;
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
	private File publicRoot; 
	private String publicPath;
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
    		this.publicPath = configuration.getString("filebroker", "public-data-path");

    		// boot up file server
    		JettyFileServer fileServer = new JettyFileServer(urlRepository);
    		fileServer.start(fileRepository.getPath(), port);

    		// cache clean up setup
    		String userDataPath = configuration.getString("filebroker", "user-data-path");
    		userDataRoot = new File(fileRepository, userDataPath);
    		publicRoot = new File(fileRepository, publicPath);
    		
    		cleanUpTriggerLimitPercentage = configuration.getInt("filebroker", "clean-up-trigger-limit-percentage");
    		cleanUpTargetPercentage = configuration.getInt("filebroker", "clean-up-target-percentage");
    		cleanUpMinimumFileAge = configuration.getInt("filebroker", "clean-up-minimum-file-age");
    		minimumSpaceForAcceptUpload = 1024*1024*configuration.getInt("filebroker", "minimum-space-for-accept-upload");
    		

    		// disable periodic clean up for now
//    		int cutoff = 1000 * configuration.getInt("filebroker", "file-life-time");
//    		int cleanUpFrequency = 1000 * configuration.getInt("filebroker", "clean-up-frequency");
//    		int checkFrequency = 1000 * 5;
//    		Timer t = new Timer("frontend-scheduled-tasks", true);
//    		t.schedule(new FileCleanUpTimerTask(userDataRoot, cutoff), 0, cleanUpFrequency);
//    		t.schedule(new JettyCheckTimerTask(fileServer), 0, checkFrequency);

    		// initialise messaging
    		this.endpoint = new MessagingEndpoint(this);
    		MessagingTopic filebrokerTopic = endpoint.createTopic(Topics.Name.AUTHORISED_FILEBROKER_TOPIC, AccessMode.READ);
    		filebrokerTopic.setListener(this);
    		this.managerClient = new ManagerClient(endpoint); 

    		// create keep-alive thread and register shutdown hook
    		KeepAliveShutdownHandler.init(this);

			logger.info("total space: " + FileUtils.byteCountToDisplaySize(userDataRoot.getTotalSpace()));
			logger.info("usable space: " + FileUtils.byteCountToDisplaySize(userDataRoot.getUsableSpace()));
			logger.info("cache clean up will start when usable space is less than: " + FileUtils.byteCountToDisplaySize((long) ((double)userDataRoot.getTotalSpace()*(double)(100-cleanUpTriggerLimitPercentage)/100)) + " (" + (100-cleanUpTriggerLimitPercentage) + "%)");
			logger.info("cache clean target usable space is:  " + FileUtils.byteCountToDisplaySize((long) ((double)userDataRoot.getTotalSpace()*(double)(100-cleanUpTargetPercentage)/100)) + " (" + (100-cleanUpTargetPercentage) + "%)");
			logger.info("minimum required space after upload: " + FileUtils.byteCountToDisplaySize(minimumSpaceForAcceptUpload));
			logger.info("will not clean up files newer than: " + (cleanUpMinimumFileAge/3600) + "h");
    		
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
				handlePublicUrlRequest(msg);

			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_PUBLIC_FILES_REQUEST.equals(((CommandMessage)msg).getCommand())) {
				handlePublicFilesRequest(msg);

			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_DISK_SPACE_REQUEST.equals(((CommandMessage)msg).getCommand())) {
				handleSpaceRequest((CommandMessage)msg);

			} else {
				logger.error("message " + msg.getMessageID() + " not understood");
			}
			
		} catch (Exception e) {
			logger.error(e, e);
		}
	}
	

	@Deprecated
	private void handlePublicUrlRequest(ChipsterMessage msg)
			throws MalformedURLException, JMSException {
		URL url = getPublicUrl();
		UrlMessage reply = new UrlMessage(url);
		endpoint.replyToMessage(msg, reply);
		managerClient.publicUrlRequest(msg.getUsername(), url);
	}

	private void handlePublicFilesRequest(ChipsterMessage msg)
			throws MalformedURLException, JMSException {
		List<URL> files = null;
		ChipsterMessage reply;
		try {
			files = getPublicFiles();
			reply = new UrlListMessage(files);
		} catch (IOException e) {
			reply = null;
		}
		endpoint.replyToMessage(msg, reply);
		//Disabled on Chipster2 backport
		//managerClient.publicFilesRequest(msg.getUsername(), files);
	}



	private void handleSpaceRequest(CommandMessage requestMessage) throws JMSException {
		long size = Long.parseLong(requestMessage.getNamedParameter(ParameterMessage.PARAMETER_DISK_SPACE));
		logger.debug("disk space request for " + size + " bytes");
		logger.debug("usable space is: " + userDataRoot.getUsableSpace());
		
		long usableSpaceSoftLimit =  (long) ((double)userDataRoot.getTotalSpace()*(double)(100-cleanUpTriggerLimitPercentage)/100);
		long usableSpaceHardLimit = minimumSpaceForAcceptUpload;
		long cleanUpTargetLimit = (long) ((double)userDataRoot.getTotalSpace()*(double)(100-cleanUpTargetPercentage)/100);
		
		
		// deal with the weird config case of soft limit being smaller than hard limit
		if (usableSpaceSoftLimit < usableSpaceHardLimit) {
			usableSpaceSoftLimit = usableSpaceHardLimit;
		}
		
		logger.debug("usable space soft limit is: " + usableSpaceSoftLimit);
		logger.debug("usable space hard limit is: " + usableSpaceHardLimit);
		
		boolean spaceAvailable;
		
		// space available, clean up limit will not be reached
		if (userDataRoot.getUsableSpace() - size >= usableSpaceSoftLimit) {
			logger.debug("enough space available, no need to do anything");
			spaceAvailable = true;
		} 

		// space available, clean up soft limit will be reached, hard will not be reached
		else if (userDataRoot.getUsableSpace() - size >= usableSpaceHardLimit) {
			logger.info("space request: " + FileUtils.byteCountToDisplaySize(size) + " usable: " + FileUtils.byteCountToDisplaySize(userDataRoot.getUsableSpace()) + 
					", usable space soft limit: " + FileUtils.byteCountToDisplaySize(usableSpaceSoftLimit) + " (" + (100-cleanUpTriggerLimitPercentage) + 
					"%) will be reached --> scheduling clean up");
			spaceAvailable = true;
			
			final long targetUsableSpace = size + cleanUpTargetLimit;
			// schedule clean up
			new Thread(new Runnable() {
				@Override
				public void run() {
					try {
						long cleanUpBeginTime = System.currentTimeMillis();
						logger.info("cache cleanup, target usable space: " + FileUtils.byteCountToDisplaySize(targetUsableSpace) + " (" + (100-cleanUpTargetPercentage) + "%)");
						Files.makeSpaceInDirectory(userDataRoot, targetUsableSpace, cleanUpMinimumFileAge, TimeUnit.SECONDS);
						logger.info("cache cleanup took " + (System.currentTimeMillis() - cleanUpBeginTime) + " ms, usable space now " + FileUtils.byteCountToDisplaySize(userDataRoot.getUsableSpace()));
					} catch (Exception e) {
						logger.warn("exception while cleaning cache", e);
					}
				}
			}, "chipster-fileserver-cache-cleanup").start();
		} 
		
		// will run out of usable space, try to make more immediately
		else if (userDataRoot.getUsableSpace() - size > 0){
			logger.info("space request: " + FileUtils.byteCountToDisplaySize(size) + " usable: " + FileUtils.byteCountToDisplaySize(userDataRoot.getUsableSpace()) + 
					", not enough space --> clean up immediately");

			try {
				long cleanUpBeginTime = System.currentTimeMillis();
				logger.info("cache cleanup, target usable space: " + FileUtils.byteCountToDisplaySize((size + cleanUpTargetLimit)) + 
						" (" + FileUtils.byteCountToDisplaySize(size) + " + " + FileUtils.byteCountToDisplaySize(cleanUpTargetLimit) + 
						" (" + (100-cleanUpTargetPercentage) + "%)");
				Files.makeSpaceInDirectory(userDataRoot, size + cleanUpTargetLimit, cleanUpMinimumFileAge, TimeUnit.SECONDS);
				logger.info("cache cleanup took " + (System.currentTimeMillis() - cleanUpBeginTime) + " ms, usable space now " + FileUtils.byteCountToDisplaySize(userDataRoot.getUsableSpace()));
			} catch (Exception e) {
				logger.warn("exception while cleaning cache", e);
			}
			logger.info("not accepting upload if less than " + FileUtils.byteCountToDisplaySize(minimumSpaceForAcceptUpload) + " usable space after upload");

			// check if cleaned up enough 
			if (userDataRoot.getUsableSpace() >= size + minimumSpaceForAcceptUpload ) {
				logger.info("enough space after cleaning");
				spaceAvailable = true;
			} else {
				logger.info("not enough space after cleaning");
				spaceAvailable = false;
			}
		} 
		
		// request more than total, no can do
		else {
			logger.info("space request: " + FileUtils.byteCountToDisplaySize(size) + " usable: " + FileUtils.byteCountToDisplaySize(userDataRoot.getUsableSpace()) + 
			", maximum space: " + FileUtils.byteCountToDisplaySize(userDataRoot.getTotalSpace()) + 
					", minimum usable: " + FileUtils.byteCountToDisplaySize(minimumSpaceForAcceptUpload) + 
					" --> not possible to make enough space");

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

	@Deprecated
	public URL getPublicUrl() throws MalformedURLException {
		return new URL(host + ":" + port + "/" + publicPath);		
	}
	
	public List<URL> getPublicFiles() throws IOException {
		
		List<URL> urlList = new LinkedList<URL>();
		
		addFilesRecursively(urlList, publicRoot);
				
		return urlList;
	}
	
	private void addFilesRecursively(List<URL> files, File path) throws IOException {
		
		for (File file : path.listFiles()) {
			
			if (file.isDirectory()) {
				addFilesRecursively(files, file);
				
			} else {

				String localPath = file.toURI().toString();//convert spaces to %20 etc.
				String publicRootString = publicRoot.toURI().toString();//convert spaces to %20 etc.

				String urlString = localPath.replace(publicRootString, host + ":" + port + "/" + publicPath + "/");

				files.add(new URL(urlString));
			}
		}		
	}
}
