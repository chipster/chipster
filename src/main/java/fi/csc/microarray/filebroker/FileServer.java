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
import java.util.zip.ZipException;

import javax.jms.JMSException;

import org.apache.activemq.DestinationDoesNotExistException;
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
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.Strings;
import fi.csc.microarray.util.SystemMonitorUtil;

public class FileServer extends NodeBase implements MessagingListener, DirectMessagingListener, ShutdownCallback {
	/**
	 * Logger for this class
	 */
	private static Logger logger;
	
	public static final String ERROR_QUOTA_EXCEEDED = "quota-exceeded";
	
	public static final String CACHE_PATH = "cache";
	public static final String STORAGE_PATH = "storage";

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
	
	private DiskCleanUp cacheCleanUp;
	
	private ExecutorService longRunningTaskExecutor = Executors.newCachedThreadPool();

	private int metadataPort;
	
	private ExampleSessionUpdater exampleSessionUpdater;

	private long defaultUserQuota; // MB
	private long quotaWarning; // percentage of the quota


	public static void main(String[] args) {
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

    		String exampleSessionPath = configuration.getString("filebroker", "example-session-path");
    		this.publicPath = configuration.getString("filebroker", "public-path");
    		this.host = configuration.getString("filebroker", "url");
    		this.port = configuration.getInt("filebroker", "port");    	
    		this.defaultUserQuota = configuration.getInt("filebroker", "default-user-quota");
    		this.quotaWarning = configuration.getInt("filebroker", "quota-warning");
    		// initialise filebroker areas
    		this.filebrokerAreas = new FileBrokerAreas(fileRepository, CACHE_PATH, STORAGE_PATH);
    		
    		// initialise url repository
    		this.urlRepository = new AuthorisedUrlRepository(host, port, CACHE_PATH, STORAGE_PATH);

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
    		
    		cacheRoot = new File(fileRepository, CACHE_PATH);
    		storageRoot = new File(fileRepository, STORAGE_PATH);
    		publicRoot = new File(fileRepository, publicPath);    		    		
    		
    		// cache clean up setup
    		int cleanUpTriggerLimitPercentage = configuration.getInt("filebroker", "clean-up-trigger-limit-percentage");
    		int cleanUpTargetPercentage = configuration.getInt("filebroker", "clean-up-target-percentage");
    		int cleanUpMinimumFileAge = configuration.getInt("filebroker", "clean-up-minimum-file-age");
    		long minimumSpaceForAcceptUpload = 1024l*1024l*configuration.getInt("filebroker", "minimum-space-for-accept-upload");    		    		    
    		
    		cacheCleanUp = new DiskCleanUp(cacheRoot, cleanUpTriggerLimitPercentage, cleanUpTargetPercentage, cleanUpMinimumFileAge, minimumSpaceForAcceptUpload);
    		
    		// boot up file server    		
    		URL hostURL = new URL(this.host);
    		JettyFileServer jettyFileServer;
    		if (externalFileServer != null) {
    			jettyFileServer = externalFileServer;
    		} else {
        		jettyFileServer = new JettyFileServer(urlRepository, metadataServer, cacheCleanUp);    			
    		}
    		jettyFileServer.start(fileRepository.getPath(), port, hostURL.getProtocol());

    		
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
    		
    		// load new zip example sessions and store as server sessions
    		if (exampleSessionPath != null && !exampleSessionPath.isEmpty()) {
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
    		}
    		
    		// create keep-alive thread and register shutdown hook
    		KeepAliveShutdownHandler.init(this);
    		
    		logger.info("minimum required space after upload: " + FileUtils.byteCountToDisplaySize(minimumSpaceForAcceptUpload));
    		logger.info("fileserver is up and running [" + ApplicationConstants.VERSION + "]");
    		logger.info("[mem: " + SystemMonitorUtil.getMemInfo() + "]");
			
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
		
	
	@SuppressWarnings("deprecation")
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
				
			} else if (msg instanceof CommandMessage && CommandMessage.COMMAND_LIST_STORAGE_USAGE_OF_SESSIONS.equals(((CommandMessage)msg).getCommand())) {
				handleListStorageUsageOfSessionsRequest(endpoint, (CommandMessage)msg, msg.getUsername());
			} else {
				logger.error("message " + msg.getMessageID() + " not understood");
			}
			
		} catch (Exception e) {
			logger.error(e, e);
		}
	}

	private void handleListStorageUsageOfSessionsRequest(
			MessagingEndpoint endpoint, CommandMessage msg, String username) throws SQLException, JMSException {
		
		CommandMessage requestMessage = (CommandMessage) msg;
		CommandMessage reply;

		List<String>[] sessions;
		sessions = metadataServer.getStorageUsageOfSessions(username);
		Long storageUsage = metadataServer.getStorageusageOfUser(username);

		reply = new CommandMessage();
		reply.addNamedParameter(ParameterMessage.PARAMETER_USERNAME_LIST, Strings.delimit(sessions[0], "\t"));
		reply.addNamedParameter(ParameterMessage.PARAMETER_SESSION_NAME_LIST, Strings.delimit(sessions[1], "\t"));
		reply.addNamedParameter(ParameterMessage.PARAMETER_SIZE_LIST, Strings.delimit(sessions[2], "\t"));
		reply.addNamedParameter(ParameterMessage.PARAMETER_DATE_LIST, Strings.delimit(sessions[3], "\t"));
		reply.addNamedParameter(ParameterMessage.PARAMETER_SESSION_UUID_LIST, Strings.delimit(sessions[4], "\t"));
		
		long quota;
		if (defaultUserQuota > 0) {
			quota = defaultUserQuota * 1024 * 1024;
		} else {
			quota = storageRoot.getUsableSpace();
		}
		
		long quotaWarningBytes = (long) (quota * quotaWarning / 100.0);
		
		reply.addNamedParameter(ParameterMessage.PARAMETER_QUOTA, "" + quota);
		reply.addNamedParameter(ParameterMessage.PARAMETER_QUOTA_WARNING, "" + quotaWarningBytes);
		reply.addNamedParameter(ParameterMessage.PARAMETER_SIZE, "" + storageUsage);
		
		jmsEndpoint.replyToMessage(requestMessage, reply);
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
		if (!AuthorisedUrlRepository.checkFilenameSyntax(fileId)) {
			reply = new CommandMessage(CommandMessage.COMMAND_FILE_OPERATION_DENIED);
		} else if (area == FileBrokerArea.STORAGE && !checkQuota(username, space)) {
			reply = new SuccessMessage(false, ERROR_QUOTA_EXCEEDED);
		} else {
			URL url = urlRepository.createAuthorisedUrl(fileId, useCompression, area, space);
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

	private void handleIsAvailable(MessagingEndpoint endpoint, ChipsterMessage msg) throws JMSException, SQLException, IOException {
		
		// parse request
		CommandMessage requestMessage = (CommandMessage) msg;
		String fileId = requestMessage.getNamedParameter(ParameterMessage.PARAMETER_FILE_ID);
		Long size = Long.parseLong(requestMessage.getNamedParameter(ParameterMessage.PARAMETER_SIZE));
		String checksum = requestMessage.getNamedParameter(ParameterMessage.PARAMETER_CHECKSUM);
		FileBrokerArea area = FileBrokerArea.valueOf(requestMessage.getNamedParameter(ParameterMessage.PARAMETER_AREA));

		
		ChipsterMessage reply;
		
		// check fileId
		if (!AuthorisedUrlRepository.checkFilenameSyntax(fileId)) {
			reply = new CommandMessage(CommandMessage.COMMAND_FILE_OPERATION_DENIED);
		} else {
			try {
				if (isAvailable(fileId, size, checksum, area)) {
					reply = new BooleanMessage(true);
				} else {
					reply = new BooleanMessage(false);
				}
			} catch ( ContentLengthException | ChecksumException e) {
				logger.info("corrupted data or data id collision (" + fileId + ", " + size + ", " + checksum + ")", e);
				reply = new CommandMessage(CommandMessage.COMMAND_FILE_OPERATION_FAILED);
			} catch (ChecksumParseException e) {
				throw new IOException(e);
			} catch (IOException e) {
				throw e;
			}
		}

		// send reply
		endpoint.replyToMessage(msg, reply);
	}

	private boolean isAvailable(String fileId, Long size, String checksum, FileBrokerArea area) throws SQLException, ChecksumParseException, IOException, ContentLengthException, ChecksumException {
				
		if (!filebrokerAreas.fileExists(fileId, area)) {
			return false;
		}
		
		Long sizeOnDisk = filebrokerAreas.getSize(fileId, area);
		Long sizeInDb = null;
		if (area == FileBrokerArea.STORAGE) {
			DbFile dbFile = metadataServer.fetchFile(fileId);
			if (dbFile == null) {
				return false;
			}
			sizeInDb = dbFile.getSize();
		}
		
		Md5FileUtils.verify(size, sizeOnDisk, sizeInDb);
		
		String checksumOnDisk = filebrokerAreas.getChecksum(fileId, area);
		
		Md5FileUtils.verify(checksum, checksumOnDisk);
		
		return true;
		
	}

	private void handleSpaceRequest(final MessagingEndpoint endpoint, final CommandMessage requestMessage) throws JMSException {
		final long size = Long.parseLong(requestMessage.getNamedParameter(ParameterMessage.PARAMETER_DISK_SPACE));
		logger.debug("disk space request for " + size + " bytes");
		logger.debug("usable space is: " + cacheRoot.getUsableSpace());
		
		longRunningTaskExecutor.execute(new Runnable() {
		
			@Override
			public void run() {
				try {
				// Schedule clean up if necessary. If the space request can't be 
				// satisfied immediately, this will wait until the clean up is done,
				// which may take several minutes.
				boolean spaceAvailable = cacheCleanUp.spaceRequest(size, true, null);

				// send reply
				BooleanMessage reply = new BooleanMessage(spaceAvailable);			
					endpoint.replyToMessage(requestMessage, reply);
				} catch (DestinationDoesNotExistException e) {
					// don't print stack trace
					logger.error("could not reply to space request, because client gave up during the clean-up");
				} catch (JMSException e) {
					logger.error("could not reply to space request", e);
				}
			}
		});
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
		List<String> fileIds = Arrays.asList(requestMessage.getNamedParameterAsArray(ParameterMessage.PARAMETER_FILE_ID_LIST));
		
		ChipsterMessage reply; 
		try {
									
			storeSession(username, name, sessionId, fileIds);			
			
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
			removeSession(previousSessionUuid);
		}					
	}

	private void handleRemoveSessionRequest(MessagingEndpoint endpoint, final CommandMessage requestMessage) throws JMSException {

		// parse request
		String sessionId = requestMessage.getNamedParameter(ParameterMessage.PARAMETER_SESSION_UUID);
		
		SuccessMessage reply;
		try {
			// if no uuid, try to get url, which was the old way
			if (sessionId == null) {
				URL url = new URL(requestMessage.getNamedParameter(ParameterMessage.PARAMETER_SESSION_URL));
				sessionId = IOUtils.getFilenameWithoutPath(url);
			}
					
			if (metadataServer.isUsernameAllowedToRemoveSession(requestMessage.getUsername(), sessionId)) {

				removeSession(sessionId);			
				reply = new SuccessMessage(true);
				
			} else {
				
				reply = new SuccessMessage(false, "user it not allowed to remove the session");
			}
			
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
			File dataFile = new File(storageRoot, removedFile);
			dataFile.delete();
			Md5FileUtils.removeMd5(dataFile);
		}		
	}

	private void handleMoveFromCacheToStorageRequest(final MessagingEndpoint endpoint, final CommandMessage requestMessage) throws JMSException, MalformedURLException {

		final String fileId = requestMessage.getNamedParameter(ParameterMessage.PARAMETER_FILE_ID);
		logger.debug("move request for: " + fileId);
		
		// check id and if in cache
		if (!(AuthorisedUrlRepository.checkFilenameSyntax(fileId) && filebrokerAreas.fileExists(fileId, FileBrokerArea.CACHE))) {
			endpoint.replyToMessage(requestMessage, new SuccessMessage(false));
		} else {

			// move
			longRunningTaskExecutor.execute(new Runnable() {

				@Override
				public void run() {

					ChipsterMessage reply = null;
					try {

						// check quota here also 
						if (!checkQuota(requestMessage.getUsername(), filebrokerAreas.getSize(fileId, FileBrokerArea.CACHE))) {
							reply = new SuccessMessage(false, ERROR_QUOTA_EXCEEDED);
						} else {

							// move the file
							boolean moveSuccess = filebrokerAreas.moveFromCacheToStorage(fileId);
							
							// add to db
							long size = filebrokerAreas.getSize(fileId, FileBrokerArea.STORAGE);
							metadataServer.addFile(fileId, size);
								
							if (moveSuccess) {
								reply = new SuccessMessage(true);
							} else {
								reply = new SuccessMessage(false);
								logger.warn("could not move from cache to storage: " + fileId);
							}
						}

					} catch (Exception e) {
						reply = new SuccessMessage(false);
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
	
	
	private boolean checkQuota(String username, long additionalBytes) throws SQLException {
		if (defaultUserQuota == -1) {
			logger.debug("quota limit disabled");
			return true;			
		} else { 
			Long usage = metadataServer.getStorageusageOfUser(username);
			boolean isAllowed = usage + additionalBytes <= defaultUserQuota*1024*1024;
			logger.debug("quota check passed: " + isAllowed);
			return isAllowed;
		}
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
					users = metadataServer.getStorageusageOfUsers();

					reply = new CommandMessage();
					reply.addNamedParameter(ParameterMessage.PARAMETER_USERNAME_LIST, Strings.delimit(users[0], "\t"));
					reply.addNamedParameter(ParameterMessage.PARAMETER_SIZE_LIST, Strings.delimit(users[1], "\t"));				

					jmsEndpoint.replyToMessage(requestMessage, reply);
				}

				// get sessions for user
				else if (msg instanceof CommandMessage && CommandMessage.COMMAND_LIST_STORAGE_USAGE_OF_SESSIONS.equals(((CommandMessage)msg).getCommand())) {
					String username = ((ParameterMessage)msg).getNamedParameter("username");
					handleListStorageUsageOfSessionsRequest(jmsEndpoint, (CommandMessage)msg, username);
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
				
				else if (msg instanceof CommandMessage && CommandMessage.COMMAND_GET_STATUS_REPORT.equals(((CommandMessage)msg).getCommand())) {
					
					CommandMessage requestMessage = (CommandMessage) msg;
					CommandMessage reply = new CommandMessage();
					
					
					FileServerAdminTools admin = getAdminTools();
					
					String sysStats = SystemMonitorUtil.getSystemStats(cacheRoot).systemStatsToString();
					
					String report = "";
					report += "SYSTEM\n\n";
					report += sysStats + "\n";
					report += admin.getDataBaseStatusReport() + "\n";
					report += admin.getStorageStatusReport(false) + "\n";

					reply.addNamedParameter(ParameterMessage.PARAMETER_STATUS_REPORT, report);
					
					jmsEndpoint.replyToMessage(requestMessage, reply);
				}
			} catch (Exception e) {
				logger.error(e, e);
			}
		}

		private FileServerAdminTools getAdminTools() {
			return new FileServerAdminTools(metadataServer, cacheRoot, storageRoot, false);
		}
	}
}
