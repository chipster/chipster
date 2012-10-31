package fi.csc.microarray.filebroker;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.jms.JMSException;

import org.apache.commons.io.FileUtils;
import org.apache.log4j.Logger;
import org.h2.tools.Server;

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
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.MemUtil;
import fi.csc.microarray.util.Strings;

public class FileServer extends NodeBase implements MessagingListener, ShutdownCallback {
	/**
	 * Logger for this class
	 */
	private static Logger logger;

	private MessagingEndpoint endpoint;	
	private ManagerClient managerClient;
	private AuthorisedUrlRepository urlRepository;
	private DerbyMetadataServer metadataServer;
	
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
    		this.publicPath = configuration.getString("filebroker", "public-path");
    		this.host = configuration.getString("filebroker", "url");
    		this.port = configuration.getInt("filebroker", "port");
    		
    		this.urlRepository = new AuthorisedUrlRepository(host, port, cachePath, storagePath);

    		// initialise metadata database
    		this.metadataServer = new DerbyMetadataServer();
    		Server h2WebConsoleServer;
    		h2WebConsoleServer = Server.createWebServer(new String[] {"-webAllowOthers",  "-webPort", String.valueOf(8082)});
    		h2WebConsoleServer.start();
    		logger.info(h2WebConsoleServer.getStatus());
    		
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

			logger.info("total space: " + FileUtils.byteCountToDisplaySize(cacheRoot.getTotalSpace()));
			logger.info("usable space: " + FileUtils.byteCountToDisplaySize(cacheRoot.getUsableSpace()));
			logger.info("cache clean up will start when usable space is less than: " + FileUtils.byteCountToDisplaySize((long) ((double)cacheRoot.getTotalSpace()*(double)(100-cleanUpTriggerLimitPercentage)/100)) + " (" + (100-cleanUpTriggerLimitPercentage) + "%)");
			logger.info("cache clean target usable space is:  " + FileUtils.byteCountToDisplaySize((long) ((double)cacheRoot.getTotalSpace()*(double)(100-cleanUpTargetPercentage)/100)) + " (" + (100-cleanUpTargetPercentage) + "%)");
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
				
				// parse request
				CommandMessage requestMessage = (CommandMessage) msg;
				boolean useCompression = requestMessage.getParameters().contains(ParameterMessage.PARAMETER_USE_COMPRESSION);
				FileBrokerArea area = FileBrokerArea.valueOf(requestMessage.getNamedParameter(ParameterMessage.PARAMETER_AREA));
				
				// check quota, if needed
				ChipsterMessage reply;
				if (area == FileBrokerArea.STORAGE && !checkQuota(msg.getUsername(), Long.parseLong(requestMessage.getNamedParameter(ParameterMessage.PARAMETER_DISK_SPACE)))) {
					reply = new CommandMessage(CommandMessage.COMMAND_FILE_OPERATION_DENIED);
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

			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_STORE_SESSION.equals(((CommandMessage)msg).getCommand())) {
				handleStoreSessionRequest((CommandMessage)msg);

			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_REMOVE_SESSION.equals(((CommandMessage)msg).getCommand())) {
				handleRemoveSessionRequest((CommandMessage)msg);

			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_LIST_SESSIONS.equals(((CommandMessage)msg).getCommand())) {
				handleListSessionsRequest((CommandMessage)msg);

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
			logger.info("space request: " + FileUtils.byteCountToDisplaySize(size) + " usable: " + FileUtils.byteCountToDisplaySize(cacheRoot.getUsableSpace()) + 
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
						Files.makeSpaceInDirectory(cacheRoot, targetUsableSpace, cleanUpMinimumFileAge, TimeUnit.SECONDS);
						logger.info("cache cleanup took " + (System.currentTimeMillis() - cleanUpBeginTime) + " ms, usable space now " + FileUtils.byteCountToDisplaySize(cacheRoot.getUsableSpace()));
					} catch (Exception e) {
						logger.warn("exception while cleaning cache", e);
					}
				}
			}, "chipster-fileserver-cache-cleanup").start();
		} 
		
		// will run out of usable space, try to make more immediately
		else if (cacheRoot.getUsableSpace() - size > 0){
			logger.info("space request: " + FileUtils.byteCountToDisplaySize(size) + " usable: " + FileUtils.byteCountToDisplaySize(cacheRoot.getUsableSpace()) + 
					", not enough space --> clean up immediately");

			try {
				long cleanUpBeginTime = System.currentTimeMillis();
				logger.info("cache cleanup, target usable space: " + FileUtils.byteCountToDisplaySize((size + cleanUpTargetLimit)) + 
						" (" + FileUtils.byteCountToDisplaySize(size) + " + " + FileUtils.byteCountToDisplaySize(cleanUpTargetLimit) + 
						" (" + (100-cleanUpTargetPercentage) + "%)");
				Files.makeSpaceInDirectory(cacheRoot, size + cleanUpTargetLimit, cleanUpMinimumFileAge, TimeUnit.SECONDS);
				logger.info("cache cleanup took " + (System.currentTimeMillis() - cleanUpBeginTime) + " ms, usable space now " + FileUtils.byteCountToDisplaySize(cacheRoot.getUsableSpace()));
			} catch (Exception e) {
				logger.warn("exception while cleaning cache", e);
			}
			logger.info("not accepting upload if less than " + FileUtils.byteCountToDisplaySize(minimumSpaceForAcceptUpload) + " usable space after upload");

			// check if cleaned up enough 
			if (cacheRoot.getUsableSpace() >= size + minimumSpaceForAcceptUpload ) {
				logger.info("enough space after cleaning");
				spaceAvailable = true;
			} else {
				logger.info("not enough space after cleaning");
				spaceAvailable = false;
			}
		} 
		
		// request more than total, no can do
		else {
			logger.info("space request: " + FileUtils.byteCountToDisplaySize(size) + " usable: " + FileUtils.byteCountToDisplaySize(cacheRoot.getUsableSpace()) + 
			", maximum space: " + FileUtils.byteCountToDisplaySize(cacheRoot.getTotalSpace()) + 
					", minimum usable: " + FileUtils.byteCountToDisplaySize(minimumSpaceForAcceptUpload) + 
					" --> not possible to make enough space");

			spaceAvailable = false;
		}
		
		// send reply
		BooleanMessage reply = new BooleanMessage(spaceAvailable);
		endpoint.replyToMessage(requestMessage, reply);
	}


	private void handleListSessionsRequest(final CommandMessage requestMessage) throws JMSException, MalformedURLException {
		String username = requestMessage.getUsername();
		CommandMessage reply;
		
		try {
			List<String> names = metadataServer.listSessionsInDatabase(username);
			reply = new CommandMessage();
			reply.addNamedParameter(ParameterMessage.PARAMETER_SESSION_NAME_LIST, Strings.delimit(names, "\t"));
			
		} catch (Exception e) {
			reply = new CommandMessage(CommandMessage.COMMAND_FILE_OPERATION_FAILED);
		}
		
		endpoint.replyToMessage(requestMessage, reply);
		
	}

	private void handleStoreSessionRequest(final CommandMessage requestMessage) throws JMSException, MalformedURLException {

		String username = requestMessage.getUsername();
		String name = requestMessage.getNamedParameter(ParameterMessage.PARAMETER_SESSION_NAME);
		URL url = urlRepository.createAuthorisedUrl(false, FileBrokerArea.STORAGE);
		
		ChipsterMessage reply; 
		try {
			metadataServer.addSessionToDatabase(username, name, IOUtils.getFilenameWithoutPath(url));
			reply = new UrlMessage(url);
			
		} catch (Exception e) {
			reply = new CommandMessage(CommandMessage.COMMAND_FILE_OPERATION_FAILED);
		}
		
		endpoint.replyToMessage(requestMessage, reply);
	}

	private void handleRemoveSessionRequest(final CommandMessage requestMessage) throws JMSException, MalformedURLException, SQLException {
		String uuid = requestMessage.getNamedParameter(ParameterMessage.PARAMETER_SESSION_UUID);
		metadataServer.removeSessionFromDatabase(uuid);
		CommandMessage reply = new CommandMessage(CommandMessage.COMMAND_FILE_OPERATION_SUCCESSFUL);
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
					reply = new CommandMessage(CommandMessage.COMMAND_FILE_OPERATION_FAILED); // TODO could add message from exception

				} catch (IOException e) {
					reply = new CommandMessage(CommandMessage.COMMAND_FILE_OPERATION_DENIED);
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



