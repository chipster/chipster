package fi.csc.microarray.filebroker;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.TimeUnit;
import java.util.zip.InflaterInputStream;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.messaging.BooleanMessageListener;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.ReplyMessageListener;
import fi.csc.microarray.messaging.UrlListMessageListener;
import fi.csc.microarray.messaging.UrlMessageListener;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.IOUtils.CopyProgressListener;
import fi.csc.microarray.util.Strings;
import fi.csc.microarray.util.UrlTransferUtil;

/**
 * Client interface for the file broker. Used by client and computing service or
 * anyone who needs transfer files within Chipster. 
 * 
 * Mostly used along the PayloadMessages which carry the URLs for the files.
 * 
 * The checkFile(URL url) method is only meant to be used to check if the cached
 * copy of a file is still available at the file broker. This method should not 
 * be used to make decisions of whether a file should be updated at the file broker.
 * In other words, the user of this class should first figure out if a file has been
 * locally modified and needs to be updated. If so, the addFile should be used 
 * directly. The checkFile() should be used if the local file has not been modified,
 * to figure out if a copy of the unmodified file still exists at the broker or should
 * a new copy of an unmodified file be added with addFile(). 
 * 
 * @author hupponen
 *
 */
public class JMSFileBrokerClient implements FileBrokerClient {
	
	private static final int SPACE_REQUEST_TIMEOUT = 300; // seconds
	private static final int QUICK_POLL_OPERATION_TIMEOUT = 5; // seconds 
	private static final int MOVE_FROM_CACHE_TO_STORAGE_TIMEOUT = 24; // hours 
	
	private static final Logger logger = Logger.getLogger(JMSFileBrokerClient.class);
	
	
	private MessagingTopic filebrokerTopic;	
	private boolean useChunked;
	private boolean useCompression;
	private String localFilebrokerPath;
	
	public JMSFileBrokerClient(MessagingTopic urlTopic, String localFilebrokerPath) throws JMSException {

		this.filebrokerTopic = urlTopic;
		this.localFilebrokerPath = localFilebrokerPath;
		
		// Read configs
		this.useChunked = DirectoryLayout.getInstance().getConfiguration().getBoolean("messaging", "use-chunked-http"); 
		this.useCompression = DirectoryLayout.getInstance().getConfiguration().getBoolean("messaging", "use-compression");
		
	}
	
	public JMSFileBrokerClient(MessagingTopic urlTopic) throws JMSException {
		this(urlTopic, null);
	}

	@Override
	public URL addSessionFile() throws JMSException, FileBrokerException {
		// quota checks are not needed for session files (metadata xml file)
		
		// return new url
		URL url = getNewUrl(useCompression, FileBrokerArea.STORAGE, 1024*1024); // assume 1 MB is enough for all session files
		if (url == null) {
			throw new FileBrokerException("filebroker is not responding");
		}
		
		return url;
	}

	/**
	 * Add file to file broker. Must be a cached file, for other types, use other versions of this method.
	 * 
	 * @see fi.csc.microarray.filebroker.FileBrokerClient#addFile(File, CopyProgressListener)
	 */
	@Override
	public URL addFile(FileBrokerArea area, File file, CopyProgressListener progressListener) throws FileBrokerException, JMSException, IOException {
		
		if (area != FileBrokerArea.CACHE) {
			throw new UnsupportedOperationException();
		}
		
		if (file.length() > 0 && !this.requestDiskSpace(file.length())) {
			throw new NotEnoughDiskSpaceException();
		}
		
		// get new url
		URL url = getNewUrl(useCompression, FileBrokerArea.CACHE, file.length());
		if (url == null) {
			throw new FileBrokerException("filebroker is not responding");
		}

		// try to move/copy it locally, or otherwise upload the file
		if (localFilebrokerPath != null && !useCompression) {
			String filename = IOUtils.getFilenameWithoutPath(url);
			File dest = new File(localFilebrokerPath, filename);
			boolean success = file.renameTo(dest);
			if (!success) {
				IOUtils.copy(file, dest); // could not move (different partition etc.), do a local copy
			}

		} else {
			InputStream stream = new FileInputStream(file);
			try {
				UrlTransferUtil.uploadStream(url, stream, useChunked, useCompression, progressListener);
			} finally {
				IOUtils.closeIfPossible(stream);
			}
		}
		
		return url;
	}



	/**
	 * @see fi.csc.microarray.filebroker.FileBrokerClient#addFile(InputStream, CopyProgressListener)
	 */
	@Override
	public URL addFile(FileBrokerArea area, InputStream file, long contentLength, CopyProgressListener progressListener) throws FileBrokerException, JMSException, IOException {
		
		URL url;
		if (area == FileBrokerArea.CACHE) {
			if (contentLength > 0  && !this.requestDiskSpace(contentLength)) {
				throw new NotEnoughDiskSpaceException();
			}

			// Get new url
			url = getNewUrl(useCompression, FileBrokerArea.CACHE, contentLength);
			if (url == null) {
				throw new FileBrokerException("New URL is null.");
			}
			
		} else {
			// Get new url
			url = getNewUrl(useCompression, FileBrokerArea.STORAGE, contentLength);
			if (url == null) {
				throw new FileBrokerException("New URL is null.");
			}
			
		}

		// Upload the stream into a file at filebroker
		logger.debug("uploading new file: " + url);
		UrlTransferUtil.uploadStream(url, file, useChunked, useCompression, progressListener);
		logger.debug("successfully uploaded: " + url);
		
		return url;
	}

	/**
	 * Add to storage a file that currently exists in cache.
	 */
	public URL addFile(InputStream file, URL cacheURL, long contentLength) throws JMSException, IOException {
		logger.debug("moving from cache to storage: " + cacheURL);

		UrlMessageListener replyListener = new UrlMessageListener();  
		URL storageURL;
		try {
			
			// ask file broker to move it or give url to upload to
			CommandMessage storeRequestMessage = new CommandMessage(CommandMessage.COMMAND_STORE_FILE);
			storeRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_URL, cacheURL.toString());
			storeRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_DISK_SPACE, Long.toString(contentLength));			
			filebrokerTopic.sendReplyableMessage(storeRequestMessage, replyListener);
			storageURL = replyListener.waitForReply(MOVE_FROM_CACHE_TO_STORAGE_TIMEOUT, TimeUnit.HOURS);
			
			// check if it was moved, if not, then upload
			if (!UrlTransferUtil.isAccessible(storageURL)) {
				// Upload the stream into a file at filebroker
				UrlTransferUtil.uploadStream(storageURL, file, useChunked, useCompression, null);				
			}
			
		} finally {
			replyListener.cleanUp();
		}
		logger.debug("storage url is: " + storageURL);
		return storageURL; 
	}


	
	/**
	 * @see fi.csc.microarray.filebroker.FileBrokerClient#getFile(java.net.URL)
	 */
	@Override
	public InputStream getFile(URL url) throws IOException {
		InputStream payload = null;
		
		// open stream
		payload = url.openStream();

		// wait for payload to become available
		long waitStartTime = System.currentTimeMillis();
		int waitTime = 10;
		while (payload.available() < 1 && (waitStartTime + QUICK_POLL_OPERATION_TIMEOUT*1000 > System.currentTimeMillis())) {
			// sleep
			try {
				Thread.sleep(waitTime);
			} catch (InterruptedException e) {
				throw new RuntimeException(e);
			}
			waitTime = waitTime*2;
		}

		// detect compression
		if (url.toString().endsWith(".compressed")) {
			return new InflaterInputStream(payload);
		} else {
			return payload;
		}
		
	}

	
	/**
	 * @see fi.csc.microarray.filebroker.FileBrokerClient#checkFile(java.net.URL, long)
	 */
	@Override
	public boolean checkFile(URL url, long contentLength) {

		HttpURLConnection connection = null;
		try {

			connection = (HttpURLConnection) url.openConnection();
		
			// check file existence
			if (connection.getResponseCode() == HttpURLConnection.HTTP_OK) {
				// should check content length and checksum
				return true;
			}
		} catch (IOException ioe) {
			return false;	
		} finally {
			IOUtils.disconnectIfPossible(connection);
		}
		return false;
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
	 */
	private URL getNewUrl(boolean useCompression, FileBrokerArea area, long contentLength) throws JMSException {
		logger.debug("getting new url");

		UrlMessageListener replyListener = new UrlMessageListener();  
		URL url;
		try {
			CommandMessage urlRequestMessage = new CommandMessage(CommandMessage.COMMAND_URL_REQUEST);
			urlRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_AREA, area.toString());
			urlRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_DISK_SPACE, Long.toString(contentLength));

			if (useCompression) {
				urlRequestMessage.addParameter(ParameterMessage.PARAMETER_USE_COMPRESSION);
			}
			filebrokerTopic.sendReplyableMessage(urlRequestMessage, replyListener);
			url = replyListener.waitForReply(SPACE_REQUEST_TIMEOUT, TimeUnit.SECONDS);
		} finally {
			replyListener.cleanUp();
		}
		logger.debug("new url is: " + url);

		return url;
	}

	/**
	 * @see fi.csc.microarray.filebroker.FileBrokerClient#getPublicFiles()
	 */
	@Override
	public List<URL> getPublicFiles() throws JMSException {
		return fetchPublicFiles();
	}

	private List<URL> fetchPublicFiles() throws JMSException {

		UrlListMessageListener replyListener = new UrlListMessageListener();  
		List<URL> urlList;
		try {
			CommandMessage fileRequestMessage = new CommandMessage(CommandMessage.COMMAND_PUBLIC_FILES_REQUEST);
			
			filebrokerTopic.sendReplyableMessage(fileRequestMessage, replyListener);
			urlList = replyListener.waitForReply(SPACE_REQUEST_TIMEOUT, TimeUnit.SECONDS);
		} finally {
			replyListener.cleanUp();
		}

		return urlList;
	}

	/**
	 * Get a local copy of a file. If the filename part of the url matches any of the files found from 
	 * local filebroker paths (given in constructor of this class), them it is symlinked or copied locally.
	 * Otherwise the file pointed by the url is downloaded.
	 * 
  	 * @see fi.csc.microarray.filebroker.FileBrokerClient#getFile(File, URL)
	 */
	@Override
	public void getFile(File destFile, URL url) throws IOException {
		
		// Try to find the file locally and symlink/copy it
		if (localFilebrokerPath != null) {
			
			// If file in filebroker cache is compressed, it will have specific suffix and we will not match it
			File fileInFilebrokerCache = new File(localFilebrokerPath, UrlTransferUtil.parseFilename(url));
			
			if (fileInFilebrokerCache.exists()) {
				boolean linkCreated = Files.createSymbolicLink(fileInFilebrokerCache, destFile);

				if (!linkCreated) {
					IOUtils.copy(fileInFilebrokerCache, destFile); // cannot create a link, must copy
				}
			}
			
		} else {
			// Not available locally, need to download
			BufferedInputStream inputStream = null;
			BufferedOutputStream fileStream = null;
			try {
				// Download to file
				inputStream = new BufferedInputStream(getFile(url));
				fileStream = new BufferedOutputStream(new FileOutputStream(destFile));
				IOUtils.copy(inputStream, fileStream);
				
			} finally {
				IOUtils.closeIfPossible(inputStream);
				IOUtils.closeIfPossible(fileStream);
			}
		}
	}

	@Override
	public boolean requestDiskSpace(long size) throws JMSException {

		BooleanMessageListener replyListener = new BooleanMessageListener();  
		Boolean spaceAvailable;
		try {
			CommandMessage spaceRequestMessage = new CommandMessage(CommandMessage.COMMAND_DISK_SPACE_REQUEST);
			spaceRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_DISK_SPACE, String.valueOf(size));
			filebrokerTopic.sendReplyableMessage(spaceRequestMessage, replyListener);
			spaceAvailable = replyListener.waitForReply(SPACE_REQUEST_TIMEOUT, TimeUnit.SECONDS);
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
	public void saveRemoteSession(String sessionName, URL sessionURL, LinkedList<URL> dataUrls) throws JMSException {
		ReplyMessageListener replyListener = new ReplyMessageListener();  
		try {
			CommandMessage storeRequestMessage = new CommandMessage(CommandMessage.COMMAND_STORE_SESSION);
			storeRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_SESSION_NAME, sessionName);
			storeRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_SESSION_UUID, sessionURL.toExternalForm());
			storeRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_FILE_URL_LIST, Strings.delimit(dataUrls, "\t"));
			
			filebrokerTopic.sendReplyableMessage(storeRequestMessage, replyListener);
			ParameterMessage reply = replyListener.waitForReply(QUICK_POLL_OPERATION_TIMEOUT, TimeUnit.SECONDS);
			
			if (reply == null || !(reply instanceof CommandMessage) || !CommandMessage.COMMAND_FILE_OPERATION_SUCCESSFUL.equals((((CommandMessage)reply).getCommand()))) {
				throw new JMSException("failed to save session metadata remotely");
			}

			
		} finally {
			replyListener.cleanUp();
		}
	}

	@Override
	public String[][] listRemoteSessions() throws JMSException {
		ReplyMessageListener replyListener = new ReplyMessageListener();  
		
		try {
			CommandMessage listRequestMessage = new CommandMessage(CommandMessage.COMMAND_LIST_SESSIONS);
			filebrokerTopic.sendReplyableMessage(listRequestMessage, replyListener);
			ParameterMessage reply = replyListener.waitForReply(QUICK_POLL_OPERATION_TIMEOUT, TimeUnit.SECONDS);
			if (reply == null) {
				throw new RuntimeException("server failed to list sessions");
			}
			String[] names, uuids;
			String namesString = reply.getNamedParameter(ParameterMessage.PARAMETER_SESSION_NAME_LIST);
			String uuidsString = reply.getNamedParameter(ParameterMessage.PARAMETER_SESSION_UUID_LIST);
			if (namesString != null && !namesString.equals("") && uuidsString != null && !uuidsString.equals("")) {
				names = namesString.split("\t");
				uuids = uuidsString.split("\t");
				if (names.length != uuids.length) {
					names = new String[0];
					uuids = new String[0];
				}
			} else {
				names = new String[0];
				uuids = new String[0];
			}
			
			return new String[][] {names, uuids};
			
		} finally {
			replyListener.cleanUp();
		}
	}

	@Override
	public void removeRemoteSession(URL sessionURL) throws JMSException {
		ReplyMessageListener replyListener = new ReplyMessageListener();  
		
		try {
			CommandMessage removeRequestMessage = new CommandMessage(CommandMessage.COMMAND_REMOVE_SESSION);
			removeRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_SESSION_UUID, sessionURL.toExternalForm()); 
			filebrokerTopic.sendReplyableMessage(removeRequestMessage, replyListener);
			ParameterMessage reply = replyListener.waitForReply(QUICK_POLL_OPERATION_TIMEOUT, TimeUnit.SECONDS);
			
			if (reply == null || !(reply instanceof CommandMessage) || !CommandMessage.COMMAND_FILE_OPERATION_SUCCESSFUL.equals((((CommandMessage)reply).getCommand()))) {
				throw new JMSException("failed to remove session");
			}
			
		} finally {
			replyListener.cleanUp();
		}
		
	}
}
