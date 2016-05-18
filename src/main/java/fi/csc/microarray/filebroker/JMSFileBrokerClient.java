package fi.csc.microarray.filebroker;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.security.KeyManagementException;
import java.security.KeyStoreException;
import java.security.NoSuchAlgorithmException;
import java.security.cert.CertificateException;
import java.util.LinkedList;
import java.util.List;
import java.util.UUID;
import java.util.concurrent.TimeUnit;
import java.util.zip.InflaterInputStream;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.messaging.AuthCancelledException;
import fi.csc.microarray.messaging.BooleanMessageListener;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.ReplyMessageListener;
import fi.csc.microarray.messaging.SuccessMessageListener;
import fi.csc.microarray.messaging.UrlListMessageListener;
import fi.csc.microarray.messaging.UrlMessageListener;
import fi.csc.microarray.messaging.admin.StorageAdminAPI.StorageEntryMessageListener;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;
import fi.csc.microarray.messaging.message.SuccessMessage;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.IOUtils.CopyProgressListener;
import fi.csc.microarray.util.KeyAndTrustManager;
import fi.csc.microarray.util.Strings;
import fi.csc.microarray.util.UrlTransferUtil;

/**
 * Client interface for the file broker. Used by client and computing service or
 * anyone who needs transfer files within Chipster. 
 * 
 * Mostly used along the PayloadMessages which carry the URLs for the files.
 * 
 * @author hupponen
 *
 */
public class JMSFileBrokerClient implements FileBrokerClient {
	
	private static final int SPACE_REQUEST_TIMEOUT = 300; // seconds
	private static final int QUICK_POLL_OPERATION_TIMEOUT = 30; // seconds
	private static final int MOVE_FROM_CACHE_TO_STORAGE_TIMEOUT = 24; // hours 
	
	private static final Logger logger = Logger.getLogger(JMSFileBrokerClient.class);
	
	private MessagingTopic filebrokerTopic;	
	private boolean useChunked;
	private boolean useCompression;
	private File localFilebrokerCache;
	private File localFilebrokerStorage;
	private boolean useChecksums;
	private String overridingFilebrokerIp;
	
	public JMSFileBrokerClient(MessagingTopic urlTopic, String localFilebrokerPath, String overridingFilebrokerIp) throws JMSException, NoSuchAlgorithmException, CertificateException, FileNotFoundException, KeyStoreException, IOException, KeyManagementException {

		this.filebrokerTopic = urlTopic;
		
		// null by default
		this.overridingFilebrokerIp = overridingFilebrokerIp;
		
		if (localFilebrokerPath != null) {
			this.localFilebrokerCache = new File(localFilebrokerPath, FileServer.CACHE_PATH);
			this.localFilebrokerStorage = new File(localFilebrokerPath, FileServer.STORAGE_PATH);
		}
		
		// Read configs
		this.useChunked = DirectoryLayout.getInstance().getConfiguration().getBoolean("messaging", "use-chunked-http"); 
		this.useCompression = DirectoryLayout.getInstance().getConfiguration().getBoolean("messaging", "use-compression");
		this.useChecksums = DirectoryLayout.getInstance().getConfiguration().getBoolean("messaging", "use-checksums");

		// Initialise keystore in case HTTPS connections are needed
		KeyAndTrustManager.initialiseTrustStore();	
	}
	
	public JMSFileBrokerClient(MessagingTopic urlTopic) throws Exception {
		this(urlTopic, null, null);
	}

	/**
	 * Add file to file broker. Must be a cached file, for other types, use other versions of this method.
	 * 
	 * @see fi.csc.microarray.filebroker.FileBrokerClient#addFile(File, CopyProgressListener)
	 */
	@Override
	public void addFile(UUID jobId, UUID sessionId, String dataId, FileBrokerArea area, File file, CopyProgressListener progressListener, String datsetName) throws FileBrokerException, IOException {
		
		if (area != FileBrokerArea.CACHE) {
			throw new UnsupportedOperationException();
		}
		
		if (file.length() > 0 && !this.requestDiskSpace(file.length())) {
			throw new NotEnoughDiskSpaceException();
		}
		
		// get new url
		URL url = getNewURL(dataId, useCompression, FileBrokerArea.CACHE, file.length());
		if (url == null) {
			throw new FileBrokerException("filebroker is not responding");
		}

		// try to move/copy it locally, or otherwise upload the file
		if (localFilebrokerCache != null && !useCompression) {
			String filename = dataId;
			File dest = new File(localFilebrokerCache, filename);
			boolean success = file.renameTo(dest);
			if (!success) {
				IOUtils.copy(file, dest); // could not move (different partition etc.), do a local copy
			}

		} else {
			InputStream stream = new FileInputStream(file);
			try {
				String md5 = UrlTransferUtil.uploadStream(url, stream, useChunked, useCompression, useChecksums, progressListener);
				logger.debug("successfully uploaded: " + url + "\tlength: " + file.length() + "\tmd5: " + md5);
			} catch (ChecksumException e) {
				// corrupted data or data id collision
				throw new IOException(e);
			} finally {
				IOUtils.closeIfPossible(stream);
			}
		}
	}
	
	/**
	 * @return md5 String of the uploaded data, if enabled in configuration
	 * @see fi.csc.microarray.filebroker.FileBrokerClient#addFile(InputStream, CopyProgressListener)
	 */
	@Override
	public String addFile(String dataId, FileBrokerArea area, InputStream file, long contentLength, CopyProgressListener progressListener) throws FileBrokerException, IOException {
		
		URL url;
		if (area == FileBrokerArea.CACHE) {
			if (contentLength > 0  && !this.requestDiskSpace(contentLength)) {
				throw new NotEnoughDiskSpaceException();
			}

			// Get new url
			url = getNewURL(dataId, useCompression, FileBrokerArea.CACHE, contentLength);
			if (url == null) {
				throw new FileBrokerException("New URL is null.");
			}
			
		} else {
			// Get new url
			url = getNewURL(dataId, useCompression, FileBrokerArea.STORAGE, contentLength);
			if (url == null) {
				throw new FileBrokerException("New URL is null.");
			}
			
		}

		// Upload the stream into a file at filebroker
		logger.debug("uploading new file: " + url);
		String md5;
		try {
			md5 = UrlTransferUtil.uploadStream(url, file, useChunked, useCompression, useChecksums, progressListener);
		} catch (ChecksumException e) {
			// corrupted data or data id collision
			throw new IOException(e);
		}
		logger.debug("successfully uploaded: " + url + "\tlength: " + contentLength + "\tmd5: " + md5);
		return md5;		
	}
	
	@Override
	public ChecksumInputStream getInputStream(String dataId) throws IOException, FileBrokerException {
		URL url = null;
		try {
			url = getURL(dataId);
		} catch (FileBrokerException e) {
			logger.error(e);
		}
		
		if (url == null) {
			throw new FileNotFoundException("file not found or filebroker didn't respond: " + dataId);
		}
		
		InputStream payload = null;

		URLConnection connection = null;
		try {
			// make sure http cache is disabled
			connection = url.openConnection();
			connection.setUseCaches(false);
			KeyAndTrustManager.configureForChipsterCertificate(connection);			
			connection.connect();

			// open stream
			payload = connection.getInputStream();

		} catch (Exception e) {
			IOUtils.closeIfPossible(payload);
			throw e;
		}

		// detect compression
		
		InputStream stream = payload;
		if (url.toString().endsWith(".compressed")) {
			stream = new InflaterInputStream(payload);
		}
		
		return new ChecksumInputStream(stream, useChecksums, connection);			
	}

	
	/**
	 * Get a local copy of a file. If the dataId  matches any of the files found from 
	 * local filebroker paths (given in constructor of this class), then it is symlinked or copied locally.
	 * Otherwise the file pointed by the dataId is downloaded.
	 * @throws JMSException 
	 * @throws ChecksumException 
	 * 
	 * @see fi.csc.microarray.filebroker.FileBrokerClient#getFile(File, URL)
	 */
	@Override
	public void getFile(UUID sessionId, String dataId, File destFile) throws IOException, FileBrokerException, ChecksumException {
		
		// Try to find the file locally and symlink/copy it
		if (localFilebrokerCache != null && localFilebrokerStorage != null) {
			
			// If file in filebroker cache is compressed, it will have specific suffix and we will not match it
			File fileInFilebrokerCache = new File(localFilebrokerCache, dataId);
			File fileInFilebrokerStorage = new File(localFilebrokerStorage, dataId);
			
			File fileInFilebroker = null;
			
			if (fileInFilebrokerCache.exists()) {
				fileInFilebroker = fileInFilebrokerCache;
			} else if (fileInFilebrokerStorage.exists()) {
				fileInFilebroker = fileInFilebrokerStorage;
			}
							
			boolean linkCreated = Files.createSymbolicLink(fileInFilebroker, destFile);

			if (!linkCreated) {
				IOUtils.copy(fileInFilebroker, destFile); // cannot create a link, must copy
			}
			
		} else {
			// Not available locally, need to download
			ChecksumInputStream inputStream = null;
			OutputStream fileStream = null;
			try {
				// Download to file
				inputStream = getInputStream(dataId);				
				fileStream = new FileOutputStream(destFile);
				
				IOUtils.copy(new BufferedInputStream(inputStream), new BufferedOutputStream(fileStream));
				
				inputStream.verifyChecksums();
				
			} finally {
				IOUtils.closeIfPossible(inputStream);
				IOUtils.closeIfPossible(fileStream);
			}
		}
	}

	@Override
	public boolean isAvailable(String dataId, Long contentLength, String checksum, FileBrokerArea area) throws FileBrokerException {
		BooleanMessageListener replyListener = new BooleanMessageListener();  
		try {
			
			String contentLengthString = null;
			if (contentLength != null) {
				contentLengthString = contentLength.toString();
			}
			
			CommandMessage requestMessage = new CommandMessage(CommandMessage.COMMAND_IS_AVAILABLE);
			requestMessage.addNamedParameter(ParameterMessage.PARAMETER_FILE_ID, dataId);
			requestMessage.addNamedParameter(ParameterMessage.PARAMETER_SIZE, contentLengthString);
			requestMessage.addNamedParameter(ParameterMessage.PARAMETER_CHECKSUM, checksum);
			requestMessage.addNamedParameter(ParameterMessage.PARAMETER_AREA, area.toString());
			filebrokerTopic.sendReplyableMessage(requestMessage, replyListener);
			
			// wait
			Boolean success = replyListener.waitForReply(QUICK_POLL_OPERATION_TIMEOUT, TimeUnit.SECONDS); 
			
			// check how it went
			
			// timeout
			if (success == null) {
				throw new RuntimeException("timeout while waiting for the filebroker");
			} else {
				return success;
			}
		} catch (JMSException | AuthCancelledException e) {
			throw new FileBrokerException(e);
		} finally {
			replyListener.cleanUp();
		}
	}

	
	@Override
	public boolean moveFromCacheToStorage(String dataId) throws FileBrokerException, AuthCancelledException {
		logger.debug("moving from cache to storage: " + dataId);
	
		SuccessMessageListener replyListener = new SuccessMessageListener();  
		try {
			
			// ask file broker to move it
			CommandMessage moveRequestMessage = new CommandMessage(CommandMessage.COMMAND_MOVE_FROM_CACHE_TO_STORAGE);
			moveRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_FILE_ID, dataId);
			filebrokerTopic.sendReplyableMessage(moveRequestMessage, replyListener);
			
			// wait
			SuccessMessage successMessage = null;
			successMessage = replyListener.waitForReply(MOVE_FROM_CACHE_TO_STORAGE_TIMEOUT, TimeUnit.HOURS);
			
			// check how it went
			
			// timeout
			if (successMessage == null) {
				throw new RuntimeException("timeout while waiting for the filebroker");
			} else if (FileServer.ERROR_QUOTA_EXCEEDED.equals(successMessage.getErrorMessage())) {
				throw new QuotaExceededException();
			} else {
				return successMessage.success();
			}
		} catch (JMSException e) {
			throw new FileBrokerException(e);
		} finally {
			replyListener.cleanUp();
		}
	}
	
	/**
	 * @see fi.csc.microarray.filebroker.FileBrokerClient#getPublicFiles()
	 */
	@Override
	public List<URL> getPublicFiles() throws FileBrokerException {
		return fetchPublicFiles();
	}

	private List<URL> fetchPublicFiles() throws FileBrokerException {

		UrlListMessageListener replyListener = new UrlListMessageListener();  
		List<URL> urlList;
		try {
			CommandMessage fileRequestMessage = new CommandMessage(CommandMessage.COMMAND_PUBLIC_FILES_REQUEST);
			
			filebrokerTopic.sendReplyableMessage(fileRequestMessage, replyListener);
			urlList = replyListener.waitForReply(QUICK_POLL_OPERATION_TIMEOUT, TimeUnit.SECONDS);
		} catch (JMSException | AuthCancelledException e) {
			throw new FileBrokerException(e);
		} finally {
			replyListener.cleanUp();
		}

		return urlList;
	}

	@Override
	public boolean requestDiskSpace(long size) throws FileBrokerException {

		BooleanMessageListener replyListener = new BooleanMessageListener();  
		Boolean spaceAvailable;
		try {
			CommandMessage spaceRequestMessage = new CommandMessage(CommandMessage.COMMAND_DISK_SPACE_REQUEST);
			spaceRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_DISK_SPACE, String.valueOf(size));
			filebrokerTopic.sendReplyableMessage(spaceRequestMessage, replyListener);
			spaceAvailable = replyListener.waitForReply(SPACE_REQUEST_TIMEOUT, TimeUnit.SECONDS);
		} catch (JMSException | AuthCancelledException e) {
			throw new FileBrokerException(e);
		} finally {
			replyListener.cleanUp();
		}

		if (spaceAvailable == null) {
			logger.warn("did not get response for space request");
			return false;
		}
		return spaceAvailable;
		
	}

	@Override
	public void saveRemoteSession(String sessionName, String sessionId, LinkedList<String> dataIds) throws FileBrokerException, AuthCancelledException {
		ReplyMessageListener replyListener = new ReplyMessageListener();  
		try {
			CommandMessage storeRequestMessage = new CommandMessage(CommandMessage.COMMAND_STORE_SESSION);
			storeRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_SESSION_NAME, sessionName);
			storeRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_SESSION_UUID, sessionId);
			storeRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_FILE_ID_LIST, Strings.delimit(dataIds, "\t"));
			
			filebrokerTopic.sendReplyableMessage(storeRequestMessage, replyListener);
			ParameterMessage reply = replyListener.waitForReply(QUICK_POLL_OPERATION_TIMEOUT, TimeUnit.SECONDS);
			
			if (reply == null || !(reply instanceof CommandMessage) || !CommandMessage.COMMAND_FILE_OPERATION_SUCCESSFUL.equals((((CommandMessage)reply).getCommand()))) {
				throw new JMSException("failed to save session metadata remotely");
			}

			
		} catch (JMSException e) {
			throw new FileBrokerException(e);
		} finally {
			replyListener.cleanUp();
		}
	}

	@Override
	public List<DbSession> listRemoteSessions() throws FileBrokerException, AuthCancelledException {
		ReplyMessageListener replyListener = new ReplyMessageListener();  
		
		try {
			CommandMessage listRequestMessage = new CommandMessage(CommandMessage.COMMAND_LIST_SESSIONS);
			filebrokerTopic.sendReplyableMessage(listRequestMessage, replyListener);
			ParameterMessage reply = replyListener.waitForReply(QUICK_POLL_OPERATION_TIMEOUT, TimeUnit.SECONDS);
			if (reply == null) {
				throw new RuntimeException("server failed to list sessions");
			}
			String[] names, sessionIds;
			String namesString = reply.getNamedParameter(ParameterMessage.PARAMETER_SESSION_NAME_LIST);
			String sessionIdsString = reply.getNamedParameter(ParameterMessage.PARAMETER_SESSION_UUID_LIST);
			
			List<DbSession> sessions = new LinkedList<>();
			
			if (namesString != null && !namesString.equals("") && sessionIdsString != null && !sessionIdsString.equals("")) {
				
				names = namesString.split("\t");
				sessionIds = sessionIdsString.split("\t");
				
				for (int i = 0; i < names.length && i < sessionIds.length; i++) {
					sessions.add(new DbSession(sessionIds[i], names[i], null));
				}
				
				if (names.length != sessionIds.length) {
					sessions.clear();
				}
			}
			
			return sessions;
			
		} catch (JMSException e) {
			throw new FileBrokerException(e);
		} finally {
			replyListener.cleanUp();
		}
	}
	
	@Override
	public StorageEntryMessageListener getStorageUsage() throws InterruptedException, FileBrokerException {
		
		StorageEntryMessageListener listener = new StorageEntryMessageListener();
		try {
			listener.query(filebrokerTopic, null);
		} catch (JMSException | AuthCancelledException e) {
			throw new FileBrokerException(e);
		}
		return listener;
	}

	@Override
	public List<DbSession> listPublicRemoteSessions() throws FileBrokerException, AuthCancelledException {
		List<DbSession> allSessions = listRemoteSessions();
		List<DbSession> publicSessions = new LinkedList<>();
		
		for (DbSession session : allSessions) {
			if (session.getName().startsWith(DerbyMetadataServer.DEFAULT_EXAMPLE_SESSION_FOLDER)) {
				publicSessions.add(session);
			}
		}		
		return publicSessions;
	}

	@Override
	public void removeRemoteSession(String dataId) throws FileBrokerException, AuthCancelledException {
		SuccessMessageListener replyListener = new SuccessMessageListener();  
		
		try {
			CommandMessage removeRequestMessage = new CommandMessage(CommandMessage.COMMAND_REMOVE_SESSION);
			removeRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_SESSION_UUID, dataId); 
			filebrokerTopic.sendReplyableMessage(removeRequestMessage, replyListener);
			SuccessMessage reply = replyListener.waitForReply(QUICK_POLL_OPERATION_TIMEOUT, TimeUnit.SECONDS);
			
			if (reply == null || !reply.success()) {
				throw new JMSException("failed to remove session");
			}
			
		} catch (JMSException e) {
			throw new FileBrokerException(e);
		} finally {
			replyListener.cleanUp();
		}
	}
	
	/**
	 * Get an URL for a new file from the file broker.
	 * 
	 * Talks to the file broker using JMS.
	 * 
	 * If useCompression is true, request an url ending with .compressed.
	 * NOTE! Compression does not work with files larger than 4 gigabytes
	 * in JDK 1.6 and earlier.
	 *  
	 * @return the new URL, may be null if file broker sends null or
	 * if reply is not received before timeout
	 * 
	 * @throws JMSException
	 * @throws FileBrokerException 
	 * @throws MalformedURLException 
	 */
	private URL getNewURL(String dataId, boolean useCompression, FileBrokerArea area, long contentLength) throws FileBrokerException, MalformedURLException {
		logger.debug("getting new url");
	
		UrlMessageListener replyListener = new UrlMessageListener();  
		URL url;
		try {
			CommandMessage urlRequestMessage = new CommandMessage(CommandMessage.COMMAND_NEW_URL_REQUEST);
			urlRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_FILE_ID, dataId);
			urlRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_AREA, area.toString());
			urlRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_DISK_SPACE, Long.toString(contentLength));
	
			if (useCompression) {
				urlRequestMessage.addParameter(ParameterMessage.PARAMETER_USE_COMPRESSION);
			}
			filebrokerTopic.sendReplyableMessage(urlRequestMessage, replyListener);
			url = replyListener.waitForReply(SPACE_REQUEST_TIMEOUT, TimeUnit.SECONDS);
		} catch (JMSException | AuthCancelledException e) {
			throw new FileBrokerException(e);
		} finally {
			replyListener.cleanUp();
		}
		
		url = applyOverridingFilebrokerIp(url);
		
		logger.debug("new url is: " + url);
	
		return url;
	}

	private URL getURL(String dataId) throws FileBrokerException, MalformedURLException {
		
		logger.debug("getting url for dataId " + dataId);
		
		UrlMessageListener replyListener = new UrlMessageListener();  
		URL url;
		try {
			CommandMessage getURLMessage = new CommandMessage(CommandMessage.COMMAND_GET_URL);
			getURLMessage.addNamedParameter(ParameterMessage.PARAMETER_FILE_ID, dataId);
	
			filebrokerTopic.sendReplyableMessage(getURLMessage, replyListener);
			url = replyListener.waitForReply(QUICK_POLL_OPERATION_TIMEOUT, TimeUnit.SECONDS);
		} catch (JMSException | AuthCancelledException e) {
			throw new FileBrokerException(e);
		} finally {
			replyListener.cleanUp();
		}
				
		url = applyOverridingFilebrokerIp(url);
		
		logger.debug("url is: " + url);
	
		return url;
	}

	private URL applyOverridingFilebrokerIp(URL url)
			throws MalformedURLException {
		if (url != null && overridingFilebrokerIp != null) {
			logger.debug("overriding filebroker ip: " + overridingFilebrokerIp);
			url = new URL(url.getProtocol(), overridingFilebrokerIp, url.getPort(), url.getFile());						
		}		
		return url;
	}

	@Override
	public String getExternalURL(String dataId) throws FileBrokerException, MalformedURLException {
		return getURL(dataId).toExternalForm();
	}

	@Override
	public Long getContentLength(String dataId) throws IOException, FileBrokerException {
		URL url = getURL(dataId);
		if (url == null) {
			return null;
		}
		return UrlTransferUtil.getContentLength(url, true);
	}
}
