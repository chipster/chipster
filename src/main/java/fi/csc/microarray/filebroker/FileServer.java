package fi.csc.microarray.filebroker;

import java.io.File;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URL;
import java.sql.SQLException;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;
import java.util.zip.ZipException;

import javax.jms.JMSException;

import org.apache.commons.io.FileUtils;
import org.apache.log4j.Logger;
import org.h2.tools.Server;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.filebroker.FileBrokerClient.FileBrokerArea;
import fi.csc.microarray.manager.ManagerClient;
import fi.csc.microarray.messaging.DirectMessagingListener;
import fi.csc.microarray.messaging.JMSMessagingEndpoint;
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
import fi.csc.microarray.messaging.message.SuccessMessage;
import fi.csc.microarray.messaging.message.UrlListMessage;
import fi.csc.microarray.messaging.message.UrlMessage;
import fi.csc.microarray.service.KeepAliveShutdownHandler;
import fi.csc.microarray.service.ShutdownCallback;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.MemUtil;
import fi.csc.microarray.util.Strings;

public class FileServer extends NodeBase implements MessagingListener, DirectMessagingListener, ShutdownCallback {
	/**
	 * Logger for this class
	 */
	private static Logger logger;

	private MessagingEndpoint jmsEndpoint;	
	private ManagerClient managerClient;
	private AuthorisedUrlRepository urlRepository;
	private FileBrokerAreas filebrokerAreas;
	private DerbyMetadataServer metadataServer;
	
	private File cacheRoot;
	private File storageRoot;
	private File publicRoot; 
	private String publicPath;
	private String host;
	private int port;

	private int cleanUpTriggerLimitPercentage;
	private int cleanUpTargetPercentage;
	private int cleanUpMinimumFileAge;
	private long minimumSpaceForAcceptUpload;
	
	private ExecutorService longRunningTaskExecutor = Executors.newCachedThreadPool();

	private int metadataPort;
	
	private ExampleSessionUpdater exampleSessionUpdater;
	private List<FileServerListener> listeners = new LinkedList<FileServerListener>();


	public static void main(String[] args) {
		DirectoryLayout.getInstance().getConfiguration();
		new FileServer(null, null);
	}

    public FileServer(String configURL, MessagingEndpoint overriddenEndpoint) {
    	this(configURL, overriddenEndpoint, null);
    }

    public FileServer(String configURL, MessagingEndpoint overriddenEndpoint, JettyFileServer externalFileServer) {

    	try {
    		// initialise dir and logging
    		DirectoryLayout.initialiseServerLayout(
    		        Arrays.asList(new String[] {"frontend", "filebroker"}),
    		        configURL);
    		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
    		
    		logger = Logger.getLogger(FileServer.class);

    		// get configurations
    		File fileRepository = DirectoryLayout.getInstance().getFileRoot();
    		String cachePath = configuration.getString("filebroker", "cache-path");
    		String storagePath = configuration.getString("filebroker", "storage-path");
    		String exampleSessionPath = configuration.getString("filebroker", "example-session-path");
    		this.publicPath = configuration.getString("filebroker", "public-path");
    		this.host = configuration.getString("filebroker", "url");
    		this.port = configuration.getInt("filebroker", "port");    	

    		// initialise filebroker areas
    		this.filebrokerAreas = new FileBrokerAreas(fileRepository, cachePath, storagePath);
    		
    		// initialise url repository
    		this.urlRepository = new AuthorisedUrlRepository(host, port, cachePath, storagePath);

    		// initialise metadata database
    		logger.info("starting derby metadata server");
    		this.metadataPort = configuration.getInt("filebroker", "metadata-port");
    		this.metadataServer = new DerbyMetadataServer();
    		if (metadataPort > 0) {
    			Server h2WebConsoleServer;
    			h2WebConsoleServer = Server.createWebServer(new String[] {"-webAllowOthers",  "-webPort", String.valueOf(this.metadataPort)});
    			h2WebConsoleServer.start();
    			logger.info("started metadata server web interface: " + h2WebConsoleServer.getStatus());
    		} else {
    			logger.info("not starting metadata server web interface");        			
    		}
    		
    		// boot up file server    		
    		JettyFileServer jettyFileServer;
    		if (externalFileServer != null) {
    			jettyFileServer = externalFileServer;
    		} else {
        		jettyFileServer = new JettyFileServer(urlRepository, metadataServer);    			
    		}
    		jettyFileServer.start(fileRepository.getPath(), port);

    		// cache clean up setup
    		cacheRoot = new File(fileRepository, cachePath);
    		publicRoot = new File(fileRepository, publicPath);
    		
    		cleanUpTriggerLimitPercentage = configuration.getInt("filebroker", "clean-up-trigger-limit-percentage");
    		cleanUpTargetPercentage = configuration.getInt("filebroker", "clean-up-target-percentage");
    		cleanUpMinimumFileAge = configuration.getInt("filebroker", "clean-up-minimum-file-age");
    		minimumSpaceForAcceptUpload = 1024*1024*configuration.getInt("filebroker", "minimum-space-for-accept-upload");
    		
    		storageRoot = new File(fileRepository, storagePath);
    		publicRoot = new File(fileRepository, publicPath);
    		
    		// disable periodic clean up for now
//    		int cutoff = 1000 * configuration.getInt("filebroker", "file-life-time");
//    		int cleanUpFrequency = 1000 * configuration.getInt("filebroker", "clean-up-frequency");
//    		int checkFrequency = 1000 * 5;
//    		Timer t = new Timer("frontend-scheduled-tasks", true);
//    		t.schedule(new FileCleanUpTimerTask(userDataRoot, cutoff), 0, cleanUpFrequency);
//    		t.schedule(new JettyCheckTimerTask(fileServer), 0, checkFrequency);
    		    	    		
    		// initialise messaging, part 1
    		if (overriddenEndpoint != null) {
    			this.jmsEndpoint = overriddenEndpoint;
    		} else {
        		this.jmsEndpoint = new JMSMessagingEndpoint(this);
    		}
    		this.managerClient = new ManagerClient(jmsEndpoint);
    		
    		addEndpoint(jmsEndpoint);

    		MessagingTopic filebrokerAdminTopic = jmsEndpoint.createTopic(Topics.Name.FILEBROKER_ADMIN_TOPIC, AccessMode.READ);
    		filebrokerAdminTopic.setListener(new FilebrokerAdminMessageListener());
    		
    		// Load new zip example sessions and store as server sessions. 
    		// ManagerClient is not really needed in this, but must be initialized to avoid NullPointerException.
    		File exampleSessionDir = new File(fileRepository, exampleSessionPath);
    		this.exampleSessionUpdater = new ExampleSessionUpdater(this, metadataServer, exampleSessionDir);
    		
    		try {
    			this.exampleSessionUpdater.importExampleSessions();
    		} catch(ZipException e) {
    			//import failed because of the broken zip file. Probably something interrupted the session export
    			//and the situation needs to be resolved manually
    			logger.error("example session import failed", e);
    			throw e;
    		}
    		
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

	protected void addEndpoint(MessagingEndpoint endpoint) throws JMSException {
		MessagingTopic filebrokerTopic = endpoint.createTopic(Topics.Name.AUTHORISED_FILEBROKER_TOPIC, AccessMode.READ);
		filebrokerTopic.setListener(this);
	}

	public String getName() {
		return "filebroker";
	}


	public void onChipsterMessage(ChipsterMessage msg) {
		onChipsterMessage(msg, this.jmsEndpoint);
	}
		
	
	public void onChipsterMessage(ChipsterMessage msg, MessagingEndpoint endpoint) {
				
		try {

			if (msg instanceof CommandMessage && CommandMessage.COMMAND_NEW_URL_REQUEST.equals(((CommandMessage)msg).getCommand())) {				
				handleNewURLRequest(endpoint, msg);
				
			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_GET_URL.equals(((CommandMessage)msg).getCommand())) {				
				handleGetURL(endpoint, msg);
			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_IS_AVAILABLE.equals(((CommandMessage)msg).getCommand())) {				
				handleIsAvailable(endpoint, msg);
			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_PUBLIC_URL_REQUEST.equals(((CommandMessage)msg).getCommand())) {
				handlePublicUrlRequest(endpoint, msg);
				
			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_PUBLIC_FILES_REQUEST.equals(((CommandMessage)msg).getCommand())) {
				handlePublicFilesRequest(endpoint, msg);

			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_DISK_SPACE_REQUEST.equals(((CommandMessage)msg).getCommand())) {
				handleSpaceRequest(endpoint, (CommandMessage)msg);

			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_MOVE_FROM_CACHE_TO_STORAGE.equals(((CommandMessage)msg).getCommand())) {
				handleMoveFromCacheToStorageRequest(endpoint, (CommandMessage)msg);

			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_STORE_SESSION.equals(((CommandMessage)msg).getCommand())) {
				handleStoreSessionRequest(endpoint, (CommandMessage)msg);

			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_REMOVE_SESSION.equals(((CommandMessage)msg).getCommand())) {
				handleRemoveSessionRequest(endpoint, (CommandMessage)msg);

			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_LIST_SESSIONS.equals(((CommandMessage)msg).getCommand())) {
				handleListSessionsRequest(endpoint, (CommandMessage)msg);

			} else {
				logger.error("message " + msg.getMessageID() + " not understood");
			}
			
		} catch (Exception e) {
			logger.error(e, e);
		}
	}

	@Deprecated
	private void handlePublicUrlRequest(MessagingEndpoint endpoint, ChipsterMessage msg)
			throws MalformedURLException, JMSException {
		URL url = getPublicUrl();
		UrlMessage reply = new UrlMessage(url);
		endpoint.replyToMessage(msg, reply);
		managerClient.publicUrlRequest(msg.getUsername(), url);
	}
	
	private void handlePublicFilesRequest(MessagingEndpoint endpoint, ChipsterMessage msg)
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
		managerClient.publicFilesRequest(msg.getUsername(), files);
	}

	private void handleNewURLRequest(MessagingEndpoint endpoint, ChipsterMessage msg)
			throws Exception {
		
		
		// parse request
		CommandMessage requestMessage = (CommandMessage) msg;
		String fileId = requestMessage.getNamedParameter(ParameterMessage.PARAMETER_FILE_ID);
		boolean useCompression = requestMessage.getParameters().contains(ParameterMessage.PARAMETER_USE_COMPRESSION);
		FileBrokerArea area = FileBrokerArea.valueOf(requestMessage.getNamedParameter(ParameterMessage.PARAMETER_AREA));
		String username = msg.getUsername();
		long space = Long.parseLong(requestMessage.getNamedParameter(ParameterMessage.PARAMETER_DISK_SPACE));
		
		logger.debug("New url request, dataId: " + fileId);
		// check quota, if needed
		ChipsterMessage reply = createNewURLReply(fileId, username, space, useCompression, area);
		
		// send reply
		endpoint.replyToMessage(msg, reply);
	}
	
	
	
	private ChipsterMessage createNewURLReply(String fileId, String username, long space, boolean useCompression, FileBrokerArea area) throws Exception {
		ChipsterMessage reply;
		if (!AuthorisedUrlRepository.checkFilenameSyntax(fileId) || (area == FileBrokerArea.STORAGE && !checkQuota(username, space))) {
			reply = new CommandMessage(CommandMessage.COMMAND_FILE_OPERATION_DENIED);
		} else {
			URL url = urlRepository.createAuthorisedUrl(fileId, useCompression, area);
			reply = new UrlMessage(url);
			managerClient.urlRequest(username, url);
		}
		return reply;
	}

	
	private void handleGetURL(MessagingEndpoint endpoint, ChipsterMessage msg) throws MalformedURLException, JMSException {
		
		// parse request
		CommandMessage requestMessage = (CommandMessage) msg;
		String fileId = requestMessage.getNamedParameter(ParameterMessage.PARAMETER_FILE_ID);
		ChipsterMessage reply;
		URL url = null;
		
		// check fileId
		if (!AuthorisedUrlRepository.checkFilenameSyntax(fileId)) {
			reply = new CommandMessage(CommandMessage.COMMAND_FILE_OPERATION_DENIED);
		// find url
		} else if (filebrokerAreas.fileExists(fileId, FileBrokerArea.CACHE)) {
			url = urlRepository.constructCacheURL(fileId, "");
		} else if (filebrokerAreas.fileExists(fileId, FileBrokerArea.STORAGE)) {
			url = urlRepository.constructStorageURL(fileId, "");
		}

		// url may be null
		reply = new UrlMessage(url);
		
		// send reply
		endpoint.replyToMessage(msg, reply);
	}

	private void handleIsAvailable(MessagingEndpoint endpoint, ChipsterMessage msg) throws JMSException {
		
		// parse request
		CommandMessage requestMessage = (CommandMessage) msg;
		String fileId = requestMessage.getNamedParameter(ParameterMessage.PARAMETER_FILE_ID);
		FileBrokerArea area = FileBrokerArea.valueOf(requestMessage.getNamedParameter(ParameterMessage.PARAMETER_AREA));

		
		ChipsterMessage reply;
		
		// check fileId
		if (!AuthorisedUrlRepository.checkFilenameSyntax(fileId)) {
			reply = new CommandMessage(CommandMessage.COMMAND_FILE_OPERATION_DENIED);
		}

		// find file
		if (filebrokerAreas.fileExists(fileId, area)) {
			reply = new BooleanMessage(true);
		} else {
			reply = new BooleanMessage(false);
		}

		// send reply
		endpoint.replyToMessage(msg, reply);
	}

	
	
	private void handleSpaceRequest(MessagingEndpoint endpoint, CommandMessage requestMessage) throws JMSException {
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


	private void handleListSessionsRequest(MessagingEndpoint endpoint, final CommandMessage requestMessage) throws JMSException, MalformedURLException {
		String username = requestMessage.getUsername();
		CommandMessage reply;
		
		try {
			List<DbSession> sessions = metadataServer.listSessions(username);
			reply = new CommandMessage();
			
			LinkedList<String> sessionIds = new LinkedList<>();
			LinkedList<String> names = new LinkedList<>();
			
			for (DbSession session : sessions) {
				names.add(session.getName());
				sessionIds.add(session.getDataId());
			}
						
			reply.addNamedParameter(ParameterMessage.PARAMETER_SESSION_NAME_LIST, Strings.delimit(names, "\t"));
			reply.addNamedParameter(ParameterMessage.PARAMETER_SESSION_UUID_LIST, Strings.delimit(sessionIds, "\t"));
			
		} catch (Exception e) {
			reply = new CommandMessage(CommandMessage.COMMAND_FILE_OPERATION_FAILED);
		}
		
		endpoint.replyToMessage(requestMessage, reply);
	}

	private void handleStoreSessionRequest(MessagingEndpoint endpoint, final CommandMessage requestMessage) throws JMSException, MalformedURLException {
		
		String username = requestMessage.getUsername();
		String name = requestMessage.getNamedParameter(ParameterMessage.PARAMETER_SESSION_NAME);		
		String sessionId = AuthorisedUrlRepository.stripCompressionSuffix(requestMessage.getNamedParameter(ParameterMessage.PARAMETER_SESSION_UUID));
		
		ChipsterMessage reply; 
		try {
			
			List<String> fileIds = Arrays.asList(requestMessage.getNamedParameter(ParameterMessage.PARAMETER_FILE_ID_LIST).split("\t"));
			
			dispatch(new FileServerListener.BeforeStoreSession(username, name, sessionId, fileIds, endpoint));			
			storeSession(username, name, sessionId, fileIds);			
			dispatch(new FileServerListener.AfterStoreSession(username, name, sessionId, fileIds, endpoint));
			
			// everything went fine
			reply = new CommandMessage(CommandMessage.COMMAND_FILE_OPERATION_SUCCESSFUL);
			
		} catch (Exception e) {
			reply = new CommandMessage(CommandMessage.COMMAND_FILE_OPERATION_FAILED);
		}
		
		endpoint.replyToMessage(requestMessage, reply);
	}

	private void storeSession(String username, String name, String sessionId, List<String> fileIds) throws SQLException {
				
		// check if we are overwriting previous session
		String previousSessionUuid = metadataServer.fetchSession(username, name);
		if (previousSessionUuid != null) {				
			// move it aside
			metadataServer.renameSession("_" + name, previousSessionUuid);
		}
		
		// store session
		metadataServer.addSession(username, name, sessionId);
		
		// link files (they have been added when uploaded)
		for (String fileId : fileIds) {
			// check if the file is stored in this file broker
			if (filebrokerAreas.fileExists(fileId, FileBrokerArea.STORAGE)) {
				metadataServer.linkFileToSession(fileId, sessionId);
			}
		}

		// remove previous
		if (previousSessionUuid != null) {				
			metadataServer.removeSession(previousSessionUuid);
		}					
	}

	private void handleRemoveSessionRequest(MessagingEndpoint endpoint, final CommandMessage requestMessage) throws JMSException {
		
		SuccessMessage reply;
		try {

			// parse request, if no uuid, try to get url, which was the old way
			String sessionId = requestMessage.getNamedParameter(ParameterMessage.PARAMETER_SESSION_UUID);
			if (sessionId == null) {
				URL url = new URL(requestMessage.getNamedParameter(ParameterMessage.PARAMETER_SESSION_URL));
				sessionId = IOUtils.getFilenameWithoutPath(url);
			}
			
			dispatch(new FileServerListener.BeforeRemoveSession(sessionId, requestMessage.getUsername(), endpoint));			
			removeSession(sessionId);			
			dispatch(new FileServerListener.AfterRemoveSession(sessionId, requestMessage.getUsername(), endpoint));
			
			// reply
			reply = new SuccessMessage(true);
		} catch (Exception e) {
			reply = new SuccessMessage(false, e);
		}

		// send
		endpoint.replyToMessage(requestMessage, reply);
	}

	protected void removeSession(String sessionId) throws SQLException {
		
		// remove from database (including related data)
		List<String> removedFiles = metadataServer.removeSession(sessionId);

		// remove from filesystem
		for (String removedFile : removedFiles) {
			new File(storageRoot, removedFile).delete();
		}		
	}

	private void handleMoveFromCacheToStorageRequest(final MessagingEndpoint endpoint, final CommandMessage requestMessage) throws JMSException, MalformedURLException {

		final String fileId = requestMessage.getNamedParameter(ParameterMessage.PARAMETER_FILE_ID);
		logger.info("move request for: " + fileId);
		
		// check id and if in cache
		if (!(AuthorisedUrlRepository.checkFilenameSyntax(fileId) && filebrokerAreas.fileExists(fileId, FileBrokerArea.CACHE))) {
			endpoint.replyToMessage(requestMessage, new BooleanMessage(false));
		} else {

			// move
			longRunningTaskExecutor.execute(new Runnable() {

				@Override
				public void run() {

					ChipsterMessage reply = null;
					try {

						// check quota here also 

						// FIXME send quota exceeded message
						//					if (!checkQuota(requestMessage.getUsername(), cacheFile.length())) {
						//						throw new IOException("quota exceeded");
						//					}

						// move the file
						if (filebrokerAreas.moveFromCacheToStorage(fileId)) {
							reply = new BooleanMessage(true);
						} else {
							reply = new BooleanMessage(false);
							logger.warn("could not move from cache to storage: " + fileId);
						}

					} catch (Exception e) {
						reply = new BooleanMessage(false);
					}

					// send reply
					try {
						endpoint.replyToMessage(requestMessage, reply);

					} catch (JMSException e) {
						logger.error("could not send reply message", e);
					}
				}
			});
		}
	}
	
	
	private boolean checkQuota(String username, long additionalBytes) {
		return true; // always allow
	}

	public void shutdown() {
		logger.info("shutdown requested");

		// close messaging endpoint
		try {
			this.jmsEndpoint.close();
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

	private class FilebrokerAdminMessageListener implements MessagingListener {

		/* (non-Javadoc)
		 * @see fi.csc.microarray.messaging.MessagingListener#onChipsterMessage(fi.csc.microarray.messaging.message.ChipsterMessage)
		 */
		@Override
		public void onChipsterMessage(ChipsterMessage msg) {

			try {

				// get totals
				if (msg instanceof CommandMessage && CommandMessage.COMMAND_LIST_STORAGE_USAGE_OF_USERS.equals(((CommandMessage)msg).getCommand())) {

					CommandMessage requestMessage = (CommandMessage) msg;
					CommandMessage reply;

					List<String>[] users;
					users = metadataServer.getListStorageusageOfUsers();

					reply = new CommandMessage();
					reply.addNamedParameter(ParameterMessage.PARAMETER_USERNAME_LIST, Strings.delimit(users[0], "\t"));
					reply.addNamedParameter(ParameterMessage.PARAMETER_SIZE_LIST, Strings.delimit(users[1], "\t"));				

					jmsEndpoint.replyToMessage(requestMessage, reply);
				}

				// get sessions for user
				else if (msg instanceof CommandMessage && CommandMessage.COMMAND_LIST_STORAGE_USAGE_OF_SESSIONS.equals(((CommandMessage)msg).getCommand())) {
					String username = ((ParameterMessage)msg).getNamedParameter("username");
					CommandMessage requestMessage = (CommandMessage) msg;
					CommandMessage reply;

					List<String>[] sessions;
					sessions = metadataServer.listStorageUsageOfSessions(username);

					reply = new CommandMessage();
					reply.addNamedParameter(ParameterMessage.PARAMETER_USERNAME_LIST, Strings.delimit(sessions[0], "\t"));
					reply.addNamedParameter(ParameterMessage.PARAMETER_SESSION_NAME_LIST, Strings.delimit(sessions[1], "\t"));
					reply.addNamedParameter(ParameterMessage.PARAMETER_SIZE_LIST, Strings.delimit(sessions[2], "\t"));
					reply.addNamedParameter(ParameterMessage.PARAMETER_DATE_LIST, Strings.delimit(sessions[3], "\t"));
					reply.addNamedParameter(ParameterMessage.PARAMETER_SESSION_UUID_LIST, Strings.delimit(sessions[4], "\t"));
					jmsEndpoint.replyToMessage(requestMessage, reply);
				}


				// get sessions for session name
				else if (msg instanceof CommandMessage && CommandMessage.COMMAND_GET_STORAGE_USAGE_TOTALS.equals(((CommandMessage)msg).getCommand())) {
					CommandMessage requestMessage = (CommandMessage) msg;
					CommandMessage reply;

					LinkedList<String> totals = new LinkedList<String>();

					totals.add(metadataServer.getStorageUsageTotal());										
					totals.add("" + storageRoot.getUsableSpace());

					reply = new CommandMessage();
					reply.addNamedParameter(ParameterMessage.PARAMETER_SIZE_LIST, Strings.delimit(totals, "\t"));				

					jmsEndpoint.replyToMessage(requestMessage, reply);
				}

				else if (msg instanceof CommandMessage && CommandMessage.COMMAND_REMOVE_SESSION.equals(((CommandMessage)msg).getCommand())) {
					handleRemoveSessionRequest(jmsEndpoint, (CommandMessage)msg);
				}
			} catch (Exception e) {
				logger.error(e);
			}
		}
	}
	
	/**
	 * Add a local listener for FileServer events.
	 * 
	 * @param listener
	 */
	public void addListener(FileServerListener listener) {
		listeners.add(listener);
	}
	
	/**
	 * Send a FileServer event to all registered listeners.
	 * 
	 * @param event
	 */
	private void dispatch(FileServerListener.Event event) {

		for (FileServerListener listener : listeners) {
			listener.listen(event);
		}
	}
}
