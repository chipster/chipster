package fi.csc.microarray.util;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.Proxy;
import java.net.ProxySelector;
import java.net.SocketAddress;
import java.net.URI;
import java.net.URL;
import java.util.LinkedList;
import java.util.List;

import javax.jms.JMSException;

public class UrlTransferUtil {

	private static final int CHUNK_SIZE = 2048;

	public static InputStream downloadStream(URL url) throws JMSException, IOException {
		return url.openStream();		
	}

	/**
	 * Uploads a file (or similar) over HTTP.
	 *  
	 * @param url address to copy to
	 * @param fis source to copy from
	 * @param useChunked use HTTP 1.1 chunked mode?
	 * @param progressListener can be null
	 * @return
	 * @throws JMSException
	 * @throws IOException
	 */
    public static URL uploadStream(URL url, InputStream fis, boolean useChunked, IOUtils.CopyProgressListener progressListener) throws IOException {

    	HttpURLConnection connection = null;
    	try {
    		connection = (HttpURLConnection)url.openConnection(); // should use openConnection(Proxy.NO_PROXY) if it actually worked
    		connection.setRequestMethod("PUT");
    		connection.setDoOutput(true);

    		if (useChunked) {
    			// use chunked mode or otherwise URLConnection loads everything into memory
    			// (chunked mode not supported before JRE 1.5)
    			connection.setChunkedStreamingMode(CHUNK_SIZE);
    		}

    		OutputStream os = null;
    		try {
    			os = connection.getOutputStream();
    			IOUtils.copy(fis, os, progressListener);

    		} catch (IOException e) {
    			e.printStackTrace();
    			throw e;

    		} finally {
    			IOUtils.closeIfPossible(os);
    			IOUtils.closeIfPossible(fis);
    		}

    		if (!isSuccessfulCode(connection.getResponseCode())) {
    			throw new IOException("PUT was not successful: "
    					+ connection.getResponseCode() + " " + connection.getResponseMessage());
    		}
    		
    	} finally {
    		IOUtils.disconnectIfPossible(connection);
    	}

        return url;
    }
    
    private static boolean isSuccessfulCode(int responseCode) {
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
}
