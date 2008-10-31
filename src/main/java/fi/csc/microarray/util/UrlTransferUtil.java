package fi.csc.microarray.util;

import org.apache.log4j.Logger;
import org.mortbay.util.IO;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;

import javax.jms.JMSException;

public class UrlTransferUtil {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(UrlTransferUtil.class);

	private static final int CHUNK_SIZE = 2048;

	public static InputStream downloadStream(URL url) throws JMSException, IOException {
		return url.openStream();		
	}
	
    public static URL uploadStream(URL url, InputStream fis, boolean useChunked) throws JMSException, IOException {

    	HttpURLConnection connection = null;
    	try {
    		connection = (HttpURLConnection)url.openConnection();
    		connection.setRequestMethod("PUT");
    		connection.setDoOutput(true);

    		logger.debug("uploading with parameters: useChunked " + useChunked + ", CHUNK_SIZE " + CHUNK_SIZE);
    		if (useChunked) {
    			// use chunked mode or otherwise URLConnection loads everything into memory
    			// (chunked mode not supported before JRE 1.5)
    			connection.setChunkedStreamingMode(CHUNK_SIZE);
    		}

    		OutputStream os = null;
    		try {
    			os = connection.getOutputStream();
    			IO.copy(fis, os);

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
}
