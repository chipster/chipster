package fi.csc.microarray.filebroker;

import java.io.InputStream;
import java.math.BigInteger;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.HashMap;
import java.util.Random;

import javax.jms.JMSException;

import fi.csc.microarray.util.IOUtils.CopyProgressListener;

/**
 * Imitation of FileBrokerClient used for testing purposes.
 * This mock class does not require any network connections.
 * 
 * @author naktinis
 *
 */
public class FileBrokerClientMock extends JMSFileBrokerClient {
    
    // Store all files in an array of streams
    HashMap<String, InputStream> files = new HashMap<String, InputStream>();
    
    public FileBrokerClientMock() throws JMSException {
        super(null);
    }
    
    public URL addFile(InputStream content, CopyProgressListener progressListener) {
        try {
            String urlString = new BigInteger(30, new Random()).toString(32);
            files.put(urlString, content);
            return new URL("file:" + urlString);
        } catch (MalformedURLException e) {
            e.printStackTrace();
        }
        return null;
    }
    
    public InputStream getFile(URL url) {
        return files.get(url.getPath());
    }
    
    public boolean checkFile(URL url, long contentLength) {
        return true;
    }
}
