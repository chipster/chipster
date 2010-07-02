package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.TsvParser;
import fi.csc.microarray.util.IOUtils;

/**
 * Handler for files accessed directly (e.g. tab-separated files).
 * 
 * @author naktinis
 *
 */
public class ChunkDataSource extends DataSource {
    
    private RandomAccessFile raFile;
    private TsvParser fileParser;
    
    public ChunkDataSource(URL url, TsvParser fileParser) throws FileNotFoundException {
        super(url);
        this.fileParser = fileParser;
    }

    public ChunkDataSource(File file, TsvParser fileParser) throws FileNotFoundException {
        super(file);
        raFile = new RandomAccessFile(file, "r");
        this.fileParser = fileParser;
    }
    
    public ChunkDataSource(URL urlRoot, String path, TsvParser fileParser)
            throws FileNotFoundException, MalformedURLException {
        super(urlRoot, path);
        this.fileParser = fileParser;
    }

    public ChunkDataSource(File fileRoot, String path, TsvParser fileParser)
            throws FileNotFoundException, MalformedURLException {
        this(new File(fileRoot, path), fileParser);
    }
  
    public int read(long filePosition, byte[] chunk) throws IOException {       
        
        if (raFile != null) {
            raFile.seek(filePosition);
            return raFile.read(chunk);
            
        } else {
            
            HttpURLConnection connection = null;
            try {
                
                connection = (HttpURLConnection)url.openConnection();
                connection.setRequestProperty("Range", "bytes=" + filePosition + "-" + (filePosition + chunk.length));
                int bytes = connection.getInputStream().read(chunk);
                
                return bytes;
                
            } finally {
                IOUtils.disconnectIfPossible(connection);
            }
        }
        
    }

    public long length() throws IOException {
        if (raFile != null) {
            return raFile.length();
            
        } else {
            HttpURLConnection connection = null;
            try {
                connection = (HttpURLConnection)url.openConnection();
                // connection.getContentLength() returns int, which is not enough
                return Long.parseLong(connection.getHeaderField("content-length"));
            } finally {
                IOUtils.disconnectIfPossible(connection);
            }
        }
    }

    public TsvParser getFileParser() {
        return fileParser;
    }
}
