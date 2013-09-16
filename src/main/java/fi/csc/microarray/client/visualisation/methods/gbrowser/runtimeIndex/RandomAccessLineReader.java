package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.ChromosomeNormaliser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;

/**
 * Custom implementation of random access line reading, because method HttpInputStream.skip()  
 * doesn't support random access and is therefore useless with big files. File implementation
 * FileInputStream.skip() is workable, but this is used also with files because of consistency and 
 * speed. This implementation is fast for getting a few lines anywhere in the file, 
 * but sequential performance is poor over http, because requests are made sequentially at the moment.
 * 
 * Jetty server sends only about 4 kilobytes per range request.
 * Sequential throughput varies greatly with ping, from 1 MB/s in local network (1 ms ping)
 * to 10 kB/s (~50 ms ping).
 * 
 * Buffering is needed in http, but with local files it might not be necessary. Probably RandomAccessFile 
 * could be queried directly, because file reading is buffered anyway in operating system level. 
 * Currently the buffer is used also with files, but the performance penalty doesn't seem to be too bad, 
 * because this implementation is still faster than FileInputStream.skip() (but there hasn't been a direct 
 * comparison with bare RandomAccessFile).   
 * 
 * @author klemela
 */
public class RandomAccessLineReader {
	
	//Must be greater than length of longest row
	public static final int HTTP_BUFFER_SIZE = 1024*4;
	
	//File position of the buffer's first byte
	private long position = -1;
	
	private String buffer;

	//Reads random access bytes from file or http
	private ByteDataSource byteDataSource;

	public RandomAccessLineReader(DataUrl dataUrl) throws URISyntaxException, IOException {
		byteDataSource = new ByteDataSource(dataUrl);
	}
	

	public RandomAccessLineReader(URL url, BedLineParser bedLineParser,
			File outputFile, ChromosomeNormaliser chromosomeNormaliser) {
		// TODO Auto-generated constructor stub
	}


	/**
	 * Set file position (in bytes) where to start reading. Return value is false, if the
	 * requested position is outside of this file. 
	 * 
	 * @param position File position in bytes.
	 * @return False if this file doesn't contain requested location.
	 * @throws IOException
	 * @throws GBrowserException
	 */
	public boolean setPosition(long position) throws IOException, GBrowserException {
		
		//Check that position is ok
		if (position < 0 || position > length() - 1) {
			position = -1;
			return false;
		}
				
		if (buffer != null && position >= this.position && position < this.position + buffer.length()) {
			
			//The buffer is still useful
			buffer = buffer.substring((int) (position - this.position));
		} else {
			
			//The old buffer is useless
			buffer = null;
		}
		
		this.position = position;
		return true;
	}

	/**
	 * Read next line starting from the file position set with method setPosition() or 
	 * the next line after previous line returned by this method. 
	 * 
	 * The first call to this method after setPosition returns either partial line, 
	 * empty line or full line. 
	 * <li> Partial line is returned when the file position is in the middle of the line. 
	 * <li> Empty line is returned if the byte at file position is new line character and 
	 * <li> full line is returned if the file position points to first byte of the line. 
	 * 
	 *  <br><br> Return value is null when the end of file is reached.
	 *  
	 * @return
	 * @throws IOException
	 */
	public String readLine() throws IOException {
		if (buffer == null) {
			fillBuffer();
		}

		int indexOfNewLine = buffer.indexOf("\n");

		while (indexOfNewLine < 0 && position + buffer.length() < length()) {
			
			//Buffer run out
			fillBuffer();
			indexOfNewLine = buffer.indexOf("\n");
			
		}

		if (indexOfNewLine < 0) {
			return null; //End of file
		}
		
		//Get the requested line from buffer
		String line = buffer.substring(0, indexOfNewLine);
		
		//Remove requested line from buffer and mark its new position
		buffer = buffer.substring(indexOfNewLine + 1);
		position += indexOfNewLine + 1;

		return line;
	}
	
	//private int fillBufferCount = 0;
	
	/**
	 * Fill internal buffer starting from the this.position.
	 * 
	 * @throws IOException
	 */
	private void fillBuffer() throws IOException {
		
		//System.out.println("FillbufferCount: " + fillBufferCount++);
			
		byte[] bytes = new byte[HTTP_BUFFER_SIZE];
		
		if (buffer == null) {
			buffer = "";
		}
		
		long refillPosition = position + buffer.length();
		
		//The last parameter 'retry' should be disabled, because it's much more
		//efficient to retry downloading only after buffer runs out. 
		int length = byteDataSource.read(refillPosition, bytes, false);
		buffer = buffer + new String(bytes, 0, length);		
		
		//System.out.println("RandomAccessLineReader.fillBuffer() Position: " + position/1024/1024 + " MB \t Length: " + buffer.lastIndexOf("\n") + " bytes");
	}


	/**
	 * Close file.
	 * 
	 */
	public void close() {
		
		if (byteDataSource != null) {
			byteDataSource.close();
			byteDataSource = null;
		}
	}
	
	/**
	 * @return File length in bytes.
	 * @throws IOException
	 */
	public long length() throws IOException {

		return byteDataSource.length();
	}
}
