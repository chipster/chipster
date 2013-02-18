package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.ChunkDataSource;

/**
 * Custom implementation of random access line reading, because method HttpInputStream.skip()  
 * doesn't support random access and is therefore useless with big files. File implementation
 * FileInputStream.skip() is workable, but this is used also with files because of consistency and 
 * speed. This implementation is fast for getting a few lines anywhere in the file, 
 * but sequential performance is poor over http, because requests are made sequentially at the moment.
 * 
 * Jetty server sends only a few kilobytes (seems to vary from 8 to 100) per range request.
 * Sequential bandwidth varies greatly with ping, from 1 MB/s in local network (1 ms ping)
 * to 10 kB/s (~50 ms ping).
 * 
 * @author klemela
 */
public class RandomAccessLineReader implements LineReader {
	
	//Must be greater than length of longest row
	private static final int HTTP_BUFFER_SIZE = 1024*8;
	private long position = -1;
	private String buffer;

	private Long length = null;
	private ChunkDataSource chunkDataSource;
	private URL url;

	public RandomAccessLineReader(URL url) throws FileNotFoundException, URISyntaxException {
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

		String line = buffer.substring(0, indexOfNewLine);
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

		return chunkDataSource.length();
	}
}
