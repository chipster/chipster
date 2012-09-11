package fi.csc.microarray.filebroker;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.filebroker.FileBrokerClient.FileBrokerArea;
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
import fi.csc.microarray.security.CryptoKey;
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

	private File cacheRoot;
	private File storageRoot;
	private String publicPath;
	private String host;
	private int port;

	private int cleanUpTriggerLimitPercentage;
	private int cleanUpTargetPercentage;
	private int cleanUpMinimumFileAge;
	private long minimumSpaceForAcceptUpload;
	
	private ExecutorService longRunningTaskExecutor = Executors.newCachedThreadPool(); 

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
    		String cachePath = configuration.getString("filebroker", "cache-path");
    		String storagePath = configuration.getString("filebroker", "storage-path");
    		this.host = configuration.getString("filebroker", "url");
    		this.port = configuration.getInt("filebroker", "port");
    		
    		this.urlRepository = new AuthorisedUrlRepository(host, port, cachePath, storagePath);
    		this.publicPath = configuration.getString("filebroker", "public-path");

    		// boot up file server
    		JettyFileServer fileServer = new JettyFileServer(urlRepository);
    		fileServer.start(fileRepository.getPath(), port);

    		// cache clean up setup
    		cacheRoot = new File(fileRepository, cachePath);
    		cleanUpTriggerLimitPercentage = configuration.getInt("filebroker", "clean-up-trigger-limit-percentage");
    		cleanUpTargetPercentage = configuration.getInt("filebroker", "clean-up-target-percentage");
    		cleanUpMinimumFileAge = configuration.getInt("filebroker", "clean-up-minimum-file-age");
    		minimumSpaceForAcceptUpload = 1024*1024*configuration.getInt("filebroker", "minimum-space-for-accept-upload");
    		
    		storageRoot = new File(fileRepository, storagePath);
    		
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
				
				// parse request
				CommandMessage requestMessage = (CommandMessage) msg;
				boolean useCompression = requestMessage.getParameters().contains(ParameterMessage.PARAMETER_USE_COMPRESSION);
				FileBrokerArea area = FileBrokerArea.valueOf(requestMessage.getNamedParameter(ParameterMessage.PARAMETER_AREA));
				
				// check quota, if needed
				ChipsterMessage reply;
				if (area == FileBrokerArea.STORAGE && !checkQuota(msg.getUsername(), Long.parseLong(requestMessage.getNamedParameter(ParameterMessage.PARAMETER_DISK_SPACE)))) {
					reply = new CommandMessage(CommandMessage.COMMAND_MOVE_DENIED);
				} else {
					URL url = urlRepository.createAuthorisedUrl(useCompression, area);
					reply = new UrlMessage(url);
					managerClient.urlRequest(msg.getUsername(), url);
				}
				
				// send reply
				endpoint.replyToMessage(msg, reply);
				
			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_PUBLIC_URL_REQUEST.equals(((CommandMessage)msg).getCommand())) {
				URL url = getPublicUrL();
				UrlMessage reply = new UrlMessage(url);
				endpoint.replyToMessage(msg, reply);
				managerClient.publicUrlRequest(msg.getUsername(), url);

			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_DISK_SPACE_REQUEST.equals(((CommandMessage)msg).getCommand())) {
				handleSpaceRequest((CommandMessage)msg);

			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_MOVE_FROM_CACHE_TO_STORAGE.equals(((CommandMessage)msg).getCommand())) {
				handleMoveRequest((CommandMessage)msg);

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
		logger.debug("usable space is: " + cacheRoot.getUsableSpace());
		
		long usableSpaceSoftLimit =  (long) ((double)cacheRoot.getTotalSpace()*(double)(100-cleanUpTriggerLimitPercentage)/100);
		long usableSpaceHardLimit = minimumSpaceForAcceptUpload;
		long cleanUpTargetLimit = (long) ((double)cacheRoot.getTotalSpace()*(double)(100-cleanUpTargetPercentage)/100);
		
		
		// deal with the weird config case of soft limit being smaller than hard limit
		if (usableSpaceSoftLimit < usableSpaceHardLimit) {
			usableSpaceSoftLimit = usableSpaceHardLimit;
		}
		
		logger.debug("usable space soft limit is: " + usableSpaceSoftLimit);
		logger.debug("usable space hard limit is: " + usableSpaceHardLimit);
		
		boolean spaceAvailable;
		
		// space available, clean up limit will not be reached
		if (cacheRoot.getUsableSpace() - size >= usableSpaceSoftLimit) {
			logger.debug("enough space available, no need to do anything");
			spaceAvailable = true;
		} 

		// space available, clean up soft limit will be reached, hard will not be reached
		else if (cacheRoot.getUsableSpace() - size >= usableSpaceHardLimit) {
			logger.info("space available, more preferred"); 
			logger.info("requested: " + size + " usable: " + cacheRoot.getUsableSpace() + ", limit: " + usableSpaceSoftLimit);
			spaceAvailable = true;
			
			final long targetUsableSpace = size + cleanUpTargetLimit;
			// schedule clean up
			new Thread(new Runnable() {
				@Override
				public void run() {
					try {
						long cleanUpBeginTime = System.currentTimeMillis();
						logger.info("cache cleanup, target usable space: " + targetUsableSpace);
						Files.makeSpaceInDirectory(cacheRoot, targetUsableSpace, cleanUpMinimumFileAge, TimeUnit.SECONDS);
						logger.info("cache cleanup took " + (System.currentTimeMillis() - cleanUpBeginTime) + " ms");
					} catch (Exception e) {
						logger.warn("exception while cleaning cache", e);
					}
				}
			}, "chipster-fileserver-cache-cleanup").start();
		} 
		
		// hard limit will be reached, try to make more immediately
		else if (cacheRoot.getUsableSpace() - size > 0){
			logger.debug("not enough space, trying to clean");
			logger.info("requested: " + size + " usable: " + cacheRoot.getUsableSpace() + ", limit: " + usableSpaceSoftLimit);
			try {
				long cleanUpBeginTime = System.currentTimeMillis();
				logger.info("cache cleanup, target usable space: " + (size + cleanUpTargetLimit));
				Files.makeSpaceInDirectory(cacheRoot, size + cleanUpTargetLimit, cleanUpMinimumFileAge, TimeUnit.SECONDS);
				logger.info("cache cleanup took " + (System.currentTimeMillis() - cleanUpBeginTime) + " ms");
			} catch (Exception e) {
				logger.warn("exception while cleaning cache", e);
			}
			logger.info("usable after cleaning: " + cacheRoot.getUsableSpace());
			logger.info("minimum extra: " + minimumSpaceForAcceptUpload);

			// check if cleaned up enough 
			if (cacheRoot.getUsableSpace() >= size + minimumSpaceForAcceptUpload ) {
				logger.info("enough after cleaning");
				spaceAvailable = true;
			} else {
				logger.info("not enough after cleaning");
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

	
	private void handleMoveRequest(final CommandMessage requestMessage) throws JMSException, MalformedURLException {

		
		final URL cacheURL = new URL(requestMessage.getNamedParameter(ParameterMessage.PARAMETER_URL));
		logger.info("move request for: " + cacheURL);

		longRunningTaskExecutor.execute(new Runnable() {

			@Override
			public void run() {

				ChipsterMessage reply = null;
				try {
					// check that url points to our cache dir
					String[] urlPathParts = cacheURL.getPath().split("/"); 
					if (urlPathParts.length != 3 || !urlPathParts[1].equals(cacheRoot.getName()) || !CryptoKey.validateKeySyntax(urlPathParts[2])) {
						logger.info("not a valid cache url: " + cacheURL);
						throw new IllegalArgumentException("not a valid cache url: " + cacheURL);
					}

					File cacheFile = new File(cacheRoot, urlPathParts[2]);

					// check that file exists
					if (!cacheFile.exists()) {
						logger.info("cache file does not exist: " + cacheFile.getAbsolutePath());
						throw new IllegalArgumentException("cache file does not exist: " + cacheFile.getAbsolutePath());
					}

					// check quota here also
					if (!checkQuota(requestMessage.getUsername(), cacheFile.length())) {
						throw new IOException("quota exceeded");
					}

					// move the file
					String storageFileName = CryptoKey.generateRandom();
					URL storageURL = null;
					if (cacheFile.renameTo(new File(storageRoot, storageFileName))) {
						storageURL = new URL(host + ":" + port + "/" + storageRoot.getName() + "/" + storageFileName);
						reply = new UrlMessage(storageURL);

					} else {
						logger.info("could not move: " + cacheFile.getAbsolutePath() + " to " + storageURL);
						throw new IllegalArgumentException("could not move: " + cacheFile.getAbsolutePath() + " to " + storageURL);
					}

				} catch (IllegalArgumentException e) {
					reply = new CommandMessage(CommandMessage.COMMAND_MOVE_FAILED); // TODO could add message from exception

				} catch (IOException e) {
					reply = new CommandMessage(CommandMessage.COMMAND_MOVE_DENIED);
				}

				// send reply
				try {
					endpoint.replyToMessage(requestMessage, reply);
					
				} catch (JMSException e) {
					logger.error(e);
				}

			}

		});


	}
	
	
	private boolean checkQuota(String username, long additionalBytes) {
		return true;
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
		return new URL(host + ":" + port + "/" + publicPath);		
	}

}



