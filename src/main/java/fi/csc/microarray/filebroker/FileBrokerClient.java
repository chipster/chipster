package fi.csc.microarray.filebroker;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.List;

import javax.jms.JMSException;

import fi.csc.microarray.util.IOUtils.CopyProgressListener;

public interface FileBrokerClient {

	/**
	 * Ask for the filebroker to try to make certain amount of disk space available.
	 * 
	 * 
	 * @param size bytes
	 * @return true if space was available
	 */
	public boolean requestDiskSpace(long size) throws JMSException;
	
	
	/**
	 * Add InputStream as file to file broker.
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
	public abstract URL addFile(InputStream file, long contentLength, CopyProgressListener progressListener) throws NotEnoughDiskSpaceException, FileBrokerException, JMSException, IOException;

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
	public abstract URL addFile(File file, CopyProgressListener progressListener) throws NotEnoughDiskSpaceException, FileBrokerException, JMSException, IOException;

	/**
	 *  Get the InputStream for a while from the FileBroker.
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
	 * @param url
	 * @return
	 * @throws IOException
	 */
	public abstract InputStream getFile(URL url) throws IOException;

	/**
	 * Get File pointed by url to destFile. Might use local file transfer instead
	 * of downloading.
	 * 
	 * @see #getFile(URL)
	 * 
	 * @param destFile destination file that must not exist
	 * @param url source file
	 * 
	 * @throws IOException
	 */
	public abstract void getFile(File destFile, URL url) throws IOException;

	/**
	 * Check if a file exists at the file broker.
	 * 
	 * This method should only be used to check if a cached file has been
	 * removed on the server side.
	 * 
	 * This method should not be used to figuring out if a cached file should
	 * be updated or not. 
	 * 
	 * TODO Check contentLength against connection.getContentLength()
	 * TODO Checking content length may not be the best idea,
	 * especially when using compression
	 * 
	 * @param url
	 * @param contentLength
	 * @return true if file exists and contentLength matches
	 * @throws IOException 
	 */
	public abstract boolean checkFile(URL url, long contentLength);
	

	/**
	 * Retrieves the list of public files or folders from the file broker. Method blocks until result is
	 * retrieved or timeout. Talks to the file broker using JMS.
	 * 
	 * @return a list of public URLs, may be null if file broker sends null or if reply is not received before timeout
	 *  
	 * @throws JMSException
	 * @throws MalformedURLException 
	 */
	public abstract List<URL> getPublicFiles() throws JMSException, MalformedURLException;

	/**
	 * Retrieves the root of the public file area from the file broker. Method blocks until result is
	 * retrieved or timeout. Talks to the file broker using JMS.
	 * 
	 * @return the new URL, may be null if file broker sends null or if reply is not received before timeout
	 *  
	 * @throws JMSException
	 */
	public abstract URL getPublicUrl() throws Exception;

}