package fi.csc.microarray.filebroker;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.LinkedList;
import java.util.List;
import java.util.UUID;

import javax.jms.JMSException;

import fi.csc.microarray.messaging.AuthCancelledException;
import fi.csc.microarray.messaging.admin.StorageAdminAPI.StorageEntryMessageListener;
import fi.csc.microarray.util.IOUtils.CopyProgressListener;

public interface FileBrokerClient {

	public static enum FileBrokerArea {
		CACHE,
		STORAGE
	}
	

	/**
	 * Ask for the filebroker to try to make certain amount of disk space available.
	 * 
	 * 
	 * @param size bytes
	 * @return true if space was available
	 */
	public boolean requestDiskSpace(long size) throws FileBrokerException;
	
	
	/**
	 * Add InputStream as file to file broker.
	 * 
	 * @param file
	 * @param contentLength -1 if unknown
	 * @param progressListener may be null
	 * @return 
	 * @throws FileBrokerException if url from file broker is null or getting url timeouts
	 * @throws JMSException
	 * @throws IOException
	 * @throws NotEnoughDiskSpaceException
	 */
	public abstract String addFile(String dataId, FileBrokerArea area, InputStream file, long contentLength, CopyProgressListener progressListener) throws NotEnoughDiskSpaceException, FileBrokerException, IOException;

	/**
	 * Add file to file broker. Might use local transfer instead of uploading.
	 * 
	 * @see #addFile(InputStream, CopyProgressListener)
	 * 
	 * @param file
	 * @param contentLength -1 if unknown
	 * @param progressListener may be null
	 * @return url to the added file
	 * @throws FileBrokerException if url from file broker is null or getting url timeouts
	 * @throws JMSException
	 * @throws IOException
	 * @throws NotEnoughDiskSpaceException
	 */
	public abstract void addFile(UUID jobId, UUID sessionId, String dataId, FileBrokerArea area, File file, CopyProgressListener progressListener, String datsetName) throws FileBrokerException, IOException;

	/**
	 *  Get the InputStream for a file from the FileBroker.
	 *  
	 *  If payload is not available right a way, wait for a while for
	 *  the payload to become available. 
	 *  
	 *  This is useful if getFile is called immediately after addFile, 
	 *  since jetty is sometimes a bit slow to write the payload 
	 *  to the disk on the server side.
	 *  
	 *  Unfortunately, the waiting slows down getting the InputStream for
	 *  empty files. Empty files are not too common though.
	 * 
	 * @return
	 * @throws IOException
	 * @throws JMSException 
	 */
	public abstract ChecksumInputStream getInputStream(String dataId) throws IOException, FileBrokerException;
	
	/**
	 * Get File pointed by url to destFile. Might use local file transfer instead
	 * of downloading.
	 * 
	 * 
	 * @param destFile destination file that must not exist
	 * 
	 * @throws IOException
	 * @throws JMSException 
	 * @throws ChecksumException 
	 */
	public abstract void getFile(UUID sessionId, String dataId, File destFile) throws IOException, FileBrokerException, ChecksumException;	

	/**
	 * Retrieves the list of public files or folders from the file broker. Method blocks until result is
	 * retrieved or timeout. Talks to the file broker using JMS.
	 * 
	 * @return a list of public URLs, may be null if file broker sends null or if reply is not received before timeout
	 *  
	 * @throws JMSException
	 * @throws MalformedURLException 
	 */
	public abstract List<URL> getPublicFiles() throws FileBrokerException, MalformedURLException;

	/**
	 * @param name
	 * @param sessionId dataId of the session metadata file
	 * @param dataIds dataIds of other files in session
	 * @throws AuthCancelledException 
	 * @throws JMSException
	 */
	public abstract void saveRemoteSession(String name, String sessionId, LinkedList<String> dataIds) throws FileBrokerException, AuthCancelledException;
	
	/**
	 * Returns storage sessions (remote sessions) available at server. Returned array contains human readable names and corresponding URL's.
	 * First name is result[0][0] and the corresponding URL is result[1][0].
	 * 
	 * @return array of names and URL's
	 * @throws AuthCancelledException 
	 */
	public abstract List<DbSession> listRemoteSessions() throws FileBrokerException, AuthCancelledException;
	public abstract List<DbSession> listPublicRemoteSessions() throws FileBrokerException, AuthCancelledException;
	/**
	 * @param dataId dataId of the session metadata file
	 * @throws AuthCancelledException 
	 * @throws JMSException
	 */
	public void removeRemoteSession(String dataId) throws FileBrokerException, AuthCancelledException;

	/**
	 * @param dataId
	 * @param contentLength null if not available
	 * @param checksum null if not available
	 * @param area
	 * @return
	 * @throws JMSException
	 */
	public boolean isAvailable(String dataId, Long contentLength, String checksum, FileBrokerArea area) throws FileBrokerException;


	public boolean moveFromCacheToStorage(String dataId) throws FileBrokerException, AuthCancelledException;


	/**
	 * Internally client should use only dataIds instead of full URL and access data through 
	 * FileBrokerClient.getInputStream(). Use this method to only to generate links that are 
	 * needed outside Chipster and therefore don't have access to FileBrokerClient. 
	 * 
	 * @param dataId
	 * @return
	 * @throws JMSException
	 * @throws FileBrokerException
	 * @throws MalformedURLException 
	 */
	public String getExternalURL(String dataId) throws FileBrokerException, MalformedURLException;


	public Long getContentLength(String dataId) throws IOException, FileBrokerException;


	StorageEntryMessageListener getStorageUsage() throws InterruptedException, FileBrokerException;
}