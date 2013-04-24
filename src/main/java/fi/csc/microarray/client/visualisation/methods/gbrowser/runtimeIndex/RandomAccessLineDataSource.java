package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;

/**
 * Data source for reading lines from file or over http, starting from defined position. 
 * The defined position is located with random access calls, so that reading is 
 * equally fast from any part of the file, regardless of its size.
 * 
 * <p>In addition to line reading of {@link RandomAccessLineReader}, this class has a method for 
 * getting the last line of file conveniently. 
 * 
 * @author klemela
 */
public class RandomAccessLineDataSource extends DataSource {
	
	private RandomAccessLineReader lineReader;
	
	public RandomAccessLineDataSource(URL url) throws FileNotFoundException, URISyntaxException {
		super(url, null);
		
		this.lineReader = new RandomAccessLineReader(url);
	}

	@Deprecated
	public RandomAccessLineDataSource(URL url, Class<? extends AreaRequestHandler> requestHandler) throws FileNotFoundException, URISyntaxException {
		super(url, requestHandler);
		
		this.lineReader = new RandomAccessLineReader(url);
	}	
	
	/**
	 * See {@link RandomAccessLineReader#setPosition(long)} for comments.
	 * 
	 * @param position
	 * @return
	 * @throws IOException
	 * @throws GBrowserException
	 */
	public boolean setLineReaderPosition(long position) throws IOException, GBrowserException {
		
		return lineReader.setPosition(position);
	}
	
	/**
	 * See {@link RandomAccessLineReader#readLine()} for comments.
	 * 
	 * @return
	 * @throws IOException
	 */
	public String getNextLine() throws IOException {
		return lineReader.readLine();
	}
	
	/**
	 * @return Last line of the file.
	 * @throws IOException
	 * @throws GBrowserException
	 */
	public String getLastLine() throws IOException, GBrowserException {
		long maxRowLength = 1024;
		
		String line;
		
		while ((line = getLastLine(maxRowLength)) == null) {
			maxRowLength *= 2;
		}
		
		return line;
	}
	
	private String getLastLine(long maxRowLength) throws IOException, GBrowserException {

		setLineReaderPosition(length() - maxRowLength);
		
		String line;
		String lastLine = null;
		String lineBeforeLast = null;
		
		while ((line = lineReader.readLine()) != null) {
			lineBeforeLast = lastLine;
			lastLine = line;
		}
		
		if (lineBeforeLast != null) {
			return lastLine;
		} else {
			return null;
		}	
	}

	public void close() {
		lineReader.close();
	}
	
	/**
	 * See {@link RandomAccessLineReader#length()} for comments.
	 * 
	 * @return
	 * @throws IOException
	 */
	public long length() throws IOException {
		return lineReader.length();
	}
}
