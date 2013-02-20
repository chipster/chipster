package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.UnsortedDataException;

/**
 * Data source for reading lines from file or over http, starting from defined position. 
 * The defined position is located with random access calls, so that reading is 
 * equally fast from any part of the file, regardless of its size.
 * 
 * @author klemela
 *
 */
public class RandomAccessLineDataSource extends DataSource {
	
	private LineReader lineReader;

	public RandomAccessLineDataSource(URL url, Class<? extends AreaRequestHandler> requestHandler) throws FileNotFoundException, URISyntaxException {
		super(url, requestHandler);
		
		this.lineReader = new RandomAccessLineReader(url);
	}

	public RandomAccessLineDataSource(URL urlRoot, String path, Class<? extends AreaRequestHandler> requestHandler)
			throws FileNotFoundException, MalformedURLException, URISyntaxException {
		
		super(urlRoot, path, requestHandler);
		
		this.lineReader = new RandomAccessLineReader(url);
	}	
	
	public void setLineReaderPosition(long position) throws IOException, GBrowserException {
		
		lineReader.setPosition(position);
	}
	
	public String getNextLine() throws IOException {
		return lineReader.readLine();
	}
	
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
	
	public long length() throws IOException {
		return lineReader.length();
	}
	
	public static void main (String args[]) throws URISyntaxException, IOException, GBrowserException {
		
		URL fileUrl = new File(System.getProperty("user.home") + "/chipster/Homo_sapiens.GRCh37.66.gtf").toURI().toURL();
		RandomAccessLineDataSource file = new RandomAccessLineDataSource(fileUrl, GtfToFeatureConversion.class);
		
		
		URL httpUrl = new URL("http://chipster-filebroker.csc.fi:7060/public/annotations/tmp/Homo_sapiens.GRCh37.66.gtf");
		RandomAccessLineDataSource http = new RandomAccessLineDataSource(httpUrl, GtfToFeatureConversion.class);
		
		manualTest(file);
		manualTest(http);
		
		iterativeTest(file, http);
	}

	private static void iterativeTest(RandomAccessLineDataSource file,
			RandomAccessLineDataSource http) throws IOException, GBrowserException {
		
		for (int i = 0; i < 1000; i++) {
			file.setLineReaderPosition(i);
			http.setLineReaderPosition(i);
			
			testCompare(file.getNextLine(), http.getNextLine(), "" + i);
		}
		
		for (int i = 1000000; i < 1001000; i++) {
			file.setLineReaderPosition(i);
			http.setLineReaderPosition(i);
			
			testCompare(file.getNextLine(), http.getNextLine(), "" + i);
		}
		
		testCompare(file.getLastLine(), http.getLastLine(), "lastLine");
		
		System.out.println("Test done");		
	}

	private static void testCompare(String line1, String line2, String description) {
		if (!line1.equals(line2)) {
			System.out.println("Lines are different!" + description);
			System.out.println("\t" + line1);
			System.out.println("\t" + line2);
		}
	}

	private static void manualTest(RandomAccessLineDataSource file)
			throws IOException, GBrowserException {
		long t = System.currentTimeMillis();
		
		System.out.println("First 10000 lines: ");
		file.setLineReaderPosition(0);
		System.out.println("\t1: " + file.getNextLine());
		
		for (int i = 0; i < 10000; i++) {
			file.getNextLine();
		}
		System.out.println("\t10000: " + file.getNextLine());
		
		System.out.println(System.currentTimeMillis() - t + " ms ");
		t = System.currentTimeMillis();
		
		System.out.println("Seek to 100MB: ");
		file.setLineReaderPosition(1*1024*1024);
		System.out.println("\t1: " + file.getNextLine());
		

		System.out.println(System.currentTimeMillis() - t + " ms ");
		t = System.currentTimeMillis();
		
		System.out.println("Seek to 300MB: ");
		file.setLineReaderPosition(3*1024*1024);
		System.out.println("\t1: " + file.getNextLine());

		System.out.println(System.currentTimeMillis() - t + " ms ");
		t = System.currentTimeMillis();
		
		System.out.println("Last line: ");
		System.out.println("\t1: " + file.getLastLine());
		
		System.out.println(System.currentTimeMillis() - t + " ms ");
		t = System.currentTimeMillis();
	}
}
