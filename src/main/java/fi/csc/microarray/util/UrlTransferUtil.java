package fi.csc.microarray.util;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.InetAddress;
import java.net.Proxy;
import java.net.ProxySelector;
import java.net.SocketAddress;
import java.net.SocketException;
import java.net.URI;
import java.net.URL;
import java.net.URLConnection;
import java.net.UnknownHostException;
import java.util.LinkedList;
import java.util.List;
import java.util.zip.Deflater;
import java.util.zip.DeflaterOutputStream;

import javax.jms.JMSException;

import fi.csc.microarray.filebroker.ChecksumException;
import fi.csc.microarray.filebroker.ChecksumInputStream;

public class UrlTransferUtil {

	public static int HTTP_TIMEOUT_MILLISECONDS = 2000;
	private static final long POST_UPLOAD_TIMEOUT_MILLISECONDS = 500;

	private static final int CHUNK_SIZE = 2048;
	
	public static InputStream downloadStream(URL url) throws JMSException, IOException {
		URLConnection connection = url.openConnection();
		KeyAndTrustManager.configureForChipsterCertificate(connection);
		return connection.getInputStream();		
	}

	
	public static String parseFilename(URL url) {
		int start = url.getPath().contains("/") ? url.getPath().lastIndexOf("/") + 1 : url.getPath().length();
		return url.getPath().substring(start);
	}

	
	/**
	 * Uploads a file (or similar) over HTTP.
	 * NOTE! Compression does not work with files larger than 4 gigabytes
	 * in JDK 1.6 and earlier.
	 *  
	 * @param url address to copy to
	 * @param fis source to copy from
	 * @param useChunked use HTTP 1.1 chunked mode?
	 * @param progressListener can be null
	 * @return
	 * @throws JMSException
	 * @throws IOException
	 * @throws ChecksumException 
	 */
    public static String uploadStream(URL url, InputStream fis, boolean useChunked, boolean compress, boolean useChecksums, IOUtils.CopyProgressListener progressListener) throws IOException, ChecksumException {

    	HttpURLConnection connection = null;
    	String checksum = null;

    	try {
    		connection = prepareForUpload(url);

    		if (useChunked) {
    			// use chunked mode or otherwise URLConnection loads everything into memory
    			// (chunked mode not supported before JRE 1.5)
    			connection.setChunkedStreamingMode(CHUNK_SIZE);
    		}
    		
    		ChecksumInputStream is = null;
    		OutputStream os = null;
    		
    		try {
    			is = new ChecksumInputStream(fis, useChecksums, connection);    					    			
    			
    			if (compress) {
        			Deflater deflater = new Deflater(Deflater.BEST_SPEED);
        			os = new DeflaterOutputStream(connection.getOutputStream(), deflater);
    			} else {
    				os = connection.getOutputStream();	
    			}
    			
    			IOUtils.copy(is, os, progressListener);
    			
    		} catch (IOException e) {
    			e.printStackTrace();
    			throw e;

    		} finally {
    			IOUtils.closeIfPossible(is);
    			IOUtils.closeIfPossible(os);
    		}
    		
    		if (!isSuccessfulCode(connection.getResponseCode())) {
    			throw new IOException("PUT was not successful: "
    					+ connection.getResponseCode() + " " + connection.getResponseMessage());
    		}
    		
    		checksum = is.verifyChecksums();
    		
		} finally {
    		IOUtils.disconnectIfPossible(connection);
    	}

    	// Wait for upload server to make the file available, so that
    	// after this method returns we can trust the file to be accessible.
    	// Could be improved by explicitly trying to read from the URL.
    	try {
			Thread.sleep(POST_UPLOAD_TIMEOUT_MILLISECONDS);
		} catch (InterruptedException e) {
			// ignore
		}
    	
    	// May be null
    	return checksum;
    }
    
    public static boolean isSuccessfulCode(int responseCode) {
		return responseCode >= 200 && responseCode < 300; // 2xx => successful
	}
    

    /**
     * Overrides system proxy settings (JVM level) to always bypass the proxy.
     * This method must be called BEFORE any upload URL objects are created.
     * It is required because JRE does not respect at all the proxy parameter
     * given to URL.openConnection(Proxy), which would be the good solution
     * for overriding proxies for uploads.
     * 
     * @see java.net.URL#openConnection(Proxy)
     */
	public static void disableProxies() {

		ProxySelector.setDefault(new ProxySelector() {

			@Override
			public void connectFailed(URI uri, SocketAddress sa, IOException ioe) {
				// we are not interested in this
			}

			@Override
			public List<Proxy> select(URI uri) {
                LinkedList<Proxy> proxies = new LinkedList<Proxy>();
                proxies.add(Proxy.NO_PROXY);
                return proxies;
			}
	    	
	    });
	}

	public static HttpURLConnection prepareForUpload(URL url) throws IOException {
		HttpURLConnection connection = (HttpURLConnection)url.openConnection(); // should use openConnection(Proxy.NO_PROXY) if it actually worked
		KeyAndTrustManager.configureForChipsterCertificate(connection);
		connection.setRequestMethod("PUT");
		connection.setDoOutput(true);
		return connection;
	}

	public static boolean isAccessible(URL url, boolean isChipsterServer) throws IOException {
		// check the URL
		HttpURLConnection connection = (HttpURLConnection)url.openConnection();
		if (isChipsterServer) {
			KeyAndTrustManager.configureForChipsterCertificate(connection);
		} else {
			KeyAndTrustManager.configureForCACertificates(connection);
		}
		connection.setConnectTimeout(HTTP_TIMEOUT_MILLISECONDS);
		connection.connect() ; 
		return connection.getResponseCode() == HttpURLConnection.HTTP_OK;
	}


	public static Long getContentLength(URL url, boolean isChipsterServer) throws IOException {
		URLConnection connection = null;
		try {
			connection = url.openConnection();
			if (isChipsterServer) {				
				KeyAndTrustManager.configureForChipsterCertificate(connection);
			} else {
				KeyAndTrustManager.configureForCACertificates(connection);
			}
			long contentLength = connection.getContentLengthLong();

			if (contentLength >= 0) {
				return contentLength;
			} else {
				throw new IOException("content length not available: " + connection.getContent());
			}
		} finally {
			IOUtils.disconnectIfPossible(connection);
		}
	}
	
	public static boolean isLocalhost(String host) throws SocketException, UnknownHostException {
		InetAddress address = InetAddress.getByName(host);
		return address.isAnyLocalAddress() || address.isLoopbackAddress();
	}
}
