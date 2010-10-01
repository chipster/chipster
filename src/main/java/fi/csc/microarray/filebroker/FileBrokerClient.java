package fi.csc.microarray.filebroker;

import java.io.IOException;
import java.io.InputStream;
import java.net.URL;

import javax.jms.JMSException;

import fi.csc.microarray.util.IOUtils.CopyProgressListener;

public interface FileBrokerClient {

	/**
	 * Add file to file broker.
	 * 
	 * @param content
	 * @param progressListener may be null
	 * @return url to the added file
	 * @throws FileBrokerException if url from file broker is null or getting url timeouts
	 * @throws JMSException
	 * @throws IOException
	 */
	public abstract URL addFile(InputStream content, CopyProgressListener progressListener) throws FileBrokerException, JMSException,
			IOException;

	/**
	 *  Get the InputStream for a while form the FileBroker.
	 *  
	 *  If payload is not available right a way, wait for awhile for
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
	 * Retrieves the root of the public file area from the file broker. Method blocks until result is
	 * retrieved or timeout. Talks to the file broker using JMS.
	 * 
	 * @return the new URL, may be null if file broker sends null or if reply is not received before timeout
	 *  
	 * @throws JMSException
	 */
	public abstract URL getPublicUrl() throws Exception;

}