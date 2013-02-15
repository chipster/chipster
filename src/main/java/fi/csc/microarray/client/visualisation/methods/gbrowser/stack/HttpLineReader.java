package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.HttpURLConnection;
import java.net.URISyntaxException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.ChunkDataSource;
import fi.csc.microarray.util.IOUtils;

/**
 * Custom implementation of line reading over http, because HttpInputStream.skip() method 
 * doesn't support random access and is therefore useless for random access of big files. 
 * This is fast for implementation is fast getting a few lines anywhere in the file, 
 * but sequential performance is poor (maybe 10 kB/s) because Jetty server gives only a few
 * kilobytes per range request and request are made sequentially at the moment.
 * 
 * @author klemela
 */
public class HttpLineReader implements LineReader {
	
	//Must be greater than length of longest row
	private static final int HTTP_BUFFER_SIZE = 1024*2;
	private long position = -1;
	private String buffer;

	private Long length = null;
	private ChunkDataSource chunkDataSource;
	private URL url;

	public HttpLineReader(URL url) throws FileNotFoundException, URISyntaxException {
		this.url = url;
		chunkDataSource = new ChunkDataSource(url, null, null);
	}
	
	public void setPosition(long position) throws IOException {
		
		if (this.position == -1 || position < this.position || position >= this.position + buffer.length()) {
			buffer = null;
		}
		this.position = position;
	}

	public String readLine() throws IOException {
		if (buffer == null) {
			fillBuffer();
		}

		int indexOfNewLine = buffer.indexOf("\n");

		if (indexOfNewLine < 0) {
			fillBuffer();
			indexOfNewLine = buffer.indexOf("\n");
			
			if (indexOfNewLine < 0) {
				return null; //End of file
			}
		}

		String line = buffer.substring(0, indexOfNewLine - 1);
		buffer = buffer.substring(indexOfNewLine + 1);
		position += indexOfNewLine + 1;

		return line;
	}

	public void fillBuffer() throws IOException {

		byte[] bytes = new byte[HTTP_BUFFER_SIZE];
		chunkDataSource.read(position, bytes);
		buffer = new String(bytes);
	}


	public void close() {
		
		if (chunkDataSource != null) {
			chunkDataSource.close();
			chunkDataSource = null;
		}
	}
	
	public long length() throws IOException {

		if (length == null) {
			HttpURLConnection connection = null;
			try {
				connection = (HttpURLConnection)url.openConnection();
				// connection.getContentLength() returns int, which is not enough
				String string = connection.getHeaderField("content-length");
				if (string == null) {
					throw new IOException("content-length unavailable for " + url);
				}
				length = Long.parseLong(connection.getHeaderField("content-length"));
			} finally {
				IOUtils.disconnectIfPossible(connection);
			}       
		} 
		return length;
	}
}
