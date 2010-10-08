package fi.csc.microarray.filebroker;

import java.io.IOException;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.zip.InflaterInputStream;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.TempTopicMessagingListenerBase;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;
import fi.csc.microarray.messaging.message.UrlMessage;
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
		 * 
		 * 
		 * @param timeout
		 * @param unit
		 * @return
		 * @throws RuntimeException if interrupted
		 */
		public URL waitForReply(long timeout, TimeUnit unit) {
			try {
				latch.await(timeout, unit);
			} catch (InterruptedException e) {
				throw new RuntimeException(e);
			}
			return this.newUrl;
		}
	}

	
	private static final int URL_REQUEST_TIMEOUT = 30; // seconds
	private static final int FILE_AVAILABLE_TIMEOUT = 10; // seconds 
	
	private static final Logger logger = Logger.getLogger(JMSFileBrokerClient.class);
	
	
	private MessagingTopic urlTopic;	
	private boolean useChunked;
	private boolean useCompression;
	
	public JMSFileBrokerClient(MessagingTopic urlTopic) throws JMSException {
		// read configs
		this.useChunked = DirectoryLayout.getInstance().getConfiguration().getBoolean("messaging", "use-chunked-http"); 
		this.useCompression = DirectoryLayout.getInstance().getConfiguration().getBoolean("messaging", "use-compression");
		
		// initialize messaging
		this.urlTopic = urlTopic;
	}
	

	/* (non-Javadoc)
	 * @see fi.csc.microarray.filebroker.FileBrokerClient#addFile(java.io.InputStream, fi.csc.microarray.util.IOUtils.CopyProgressListener)
	 */
	public URL addFile(InputStream content, CopyProgressListener progressListener) throws FileBrokerException, JMSException, IOException {
		
		// get new url
		URL url = getNewUrl(useCompression);
		if (url == null) {
			throw new FileBrokerException("New URL is null.");
		}
		
		// upload content
		logger.debug("uploading new file: " + url);
		UrlTransferUtil.uploadStream(url, content, useChunked, useCompression, progressListener);
		logger.debug("successfully uploaded: " + url);
		
		return url;
	}
	
	
	/* (non-Javadoc)
	 * @see fi.csc.microarray.filebroker.FileBrokerClient#getFile(java.net.URL)
	 */
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

	
	/* (non-Javadoc)
	 * @see fi.csc.microarray.filebroker.FileBrokerClient#checkFile(java.net.URL, long)
	 */
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
	 * If useCompression is true, request an url ending with .compressed
	 * 
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
			url = replyListener.waitForReply(URL_REQUEST_TIMEOUT, TimeUnit.SECONDS);
		} finally {
			replyListener.cleanUp();
		}
		logger.debug("new url is: " + url);

		return url;
	}

	/* (non-Javadoc)
	 * @see fi.csc.microarray.filebroker.FileBrokerClient#getPublicUrl()
	 */
	public URL getPublicUrl() throws JMSException {

		UrlMessageListener replyListener = new UrlMessageListener();  
		URL url;
		try {
			CommandMessage urlRequestMessage = new CommandMessage(CommandMessage.COMMAND_PUBLIC_URL_REQUEST);
			urlTopic.sendReplyableMessage(urlRequestMessage, replyListener);
			url = replyListener.waitForReply(URL_REQUEST_TIMEOUT, TimeUnit.SECONDS);
		} finally {
			replyListener.cleanUp();
		}

		return url;
	}

}
