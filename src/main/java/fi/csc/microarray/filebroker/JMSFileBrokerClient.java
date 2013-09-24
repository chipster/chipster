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
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.zip.InflaterInputStream;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.messaging.BooleanMessageListener;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.TempTopicMessagingListenerBase;
import fi.csc.microarray.messaging.UrlListMessageListener;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;
import fi.csc.microarray.messaging.message.UrlMessage;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.UrlTransferUtil;
import fi.csc.microarray.util.IOUtils.CopyProgressListener;

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

	/**
	 * Reply listener for the url request.
	 * 
	 */
	private class UrlMessageListener extends TempTopicMessagingListenerBase {
		
		private URL newUrl;
		private CountDownLatch latch = new CountDownLatch(1);
		
		public void onChipsterMessage(ChipsterMessage msg) {
			if (msg instanceof UrlMessage) {
				UrlMessage urlMessage = (UrlMessage) msg;
				this.newUrl = urlMessage.getUrl();
				latch.countDown();
			}
		}
	
		/**
		 * @param timeout in given units
		 * @param unit unit of the timeout
		 * @return 
		 * @throws RuntimeException if interrupted
		 */
		public URL waitForReply(long timeout, TimeUnit unit) {
			try {
				latch.await(timeout, unit);
			} catch (InterruptedException e) {
	            logger.warn("interrupted while waiting for latch", e);
			}
			return this.newUrl;
		}
	}

	
	private static final int SPACE_REQUEST_TIMEOUT = 300; // seconds
	private static final int FILE_AVAILABLE_TIMEOUT = 5; // seconds 
	
	private static final Logger logger = Logger.getLogger(JMSFileBrokerClient.class);
	
	
	private MessagingTopic urlTopic;	
	private boolean useChunked;
	private boolean useCompression;
	private String localFilebrokerPath;
	
	public JMSFileBrokerClient(MessagingTopic urlTopic, String localFilebrokerPath) throws JMSException {

		this.urlTopic = urlTopic;
		this.localFilebrokerPath = localFilebrokerPath;
		
		// Read configs
		this.useChunked = DirectoryLayout.getInstance().getConfiguration().getBoolean("messaging", "use-chunked-http"); 
		this.useCompression = DirectoryLayout.getInstance().getConfiguration().getBoolean("messaging", "use-compression");
		
	}
	
	public JMSFileBrokerClient(MessagingTopic urlTopic) throws JMSException {
		this(urlTopic, null);
	}



	/**
	 * @see fi.csc.microarray.filebroker.FileBrokerClient#addFile(File, CopyProgressListener)
	 */
	@Override
	public URL addFile(File file, CopyProgressListener progressListener) throws FileBrokerException, JMSException, IOException {
		if (file.length() > 0 && !this.requestDiskSpace(file.length())) {
			throw new NotEnoughDiskSpaceException();
		}
		
		// Get new url
		URL url = getNewUrl(useCompression);
		if (url == null) {
			throw new FileBrokerException("New URL is null.");
		}

		// Try to move/copy it locally, or otherwise upload the file
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
	public URL addFile(InputStream file, long contentLength, CopyProgressListener progressListener) throws FileBrokerException, JMSException, IOException {
		if (contentLength > 0  && !this.requestDiskSpace(contentLength)) {
			throw new NotEnoughDiskSpaceException();
		}
		
		// Get new url
		URL url = getNewUrl(useCompression);
		if (url == null) {
			throw new FileBrokerException("New URL is null.");
		}

		// Upload the stream into a file at filebroker
		logger.debug("uploading new file: " + url);
		UrlTransferUtil.uploadStream(url, file, useChunked, useCompression, progressListener);
		logger.debug("successfully uploaded: " + url);
		
		return url;
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
		while (payload.available() < 1 && (waitStartTime + FILE_AVAILABLE_TIMEOUT*1000 > System.currentTimeMillis())) {
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
	private URL getNewUrl(boolean useCompression) throws JMSException {
		logger.debug("getting new url");

		UrlMessageListener replyListener = new UrlMessageListener();  
		URL url;
		try {
			CommandMessage urlRequestMessage = new CommandMessage(CommandMessage.COMMAND_URL_REQUEST);
			if (useCompression) {
				urlRequestMessage.addParameter(ParameterMessage.PARAMETER_USE_COMPRESSION);
			}
			urlTopic.sendReplyableMessage(urlRequestMessage, replyListener);
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
			//Chipster2 backport fix
			urlTopic.sendReplyableMessage(fileRequestMessage, replyListener);
			urlList = replyListener.waitForReply(FILE_AVAILABLE_TIMEOUT, TimeUnit.SECONDS);
		} finally {
			replyListener.cleanUp();
		}

		return urlList;
	}


	/**
	 * @see fi.csc.microarray.filebroker.FileBrokerClient#getPublicUrl()
	 */
	@Override
	public URL getPublicUrl() throws JMSException {
		return fetchPublicUrl();
	}

	private URL fetchPublicUrl() throws JMSException {

		UrlMessageListener replyListener = new UrlMessageListener();  
		URL url;
		try {
			CommandMessage urlRequestMessage = new CommandMessage(CommandMessage.COMMAND_PUBLIC_URL_REQUEST);
			urlTopic.sendReplyableMessage(urlRequestMessage, replyListener);
			url = replyListener.waitForReply(SPACE_REQUEST_TIMEOUT, TimeUnit.SECONDS);
		} finally {
			replyListener.cleanUp();
		}

		return url;
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
			urlTopic.sendReplyableMessage(spaceRequestMessage, replyListener);
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

}
