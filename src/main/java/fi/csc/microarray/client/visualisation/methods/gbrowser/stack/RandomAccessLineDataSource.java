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
 * Data source for reading lines from file or http, starting from defined position. 
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
		
		if (file != null) {
			this.lineReader = new FileLineReader(url);
		} else {
			this.lineReader = new HttpLineReader(url);
		}
	}

	public RandomAccessLineDataSource(URL urlRoot, String path, Class<? extends AreaRequestHandler> requestHandler)
			throws FileNotFoundException, MalformedURLException, URISyntaxException {
		
		super(urlRoot, path, requestHandler);
		
		if (file != null) {
			this.lineReader = new FileLineReader(url);
		} else {
			this.lineReader = new HttpLineReader(url);
		}
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

	public void checkSorting() throws IOException, UnsortedDataException {
//		byte[] bytes = new byte[100000];
//		this.read(getHeaderLength(), bytes);
//		Chunk chunk = new Chunk(new String(bytes));
//		List<RegionContent> regions = fileParser.getAll(chunk, new LinkedList<ColumnType>());
//		
//		Region previousRegion = null;
//		for (RegionContent region : regions) {
//			if (previousRegion != null) {
//				if (previousRegion.compareTo(region.region) > 0) {
//					throw new UnsortedDataException("File " + file + " isn't sorted correctly. " +
//							"Please sort the file first.");
//				}
//			}
//			previousRegion = region.region;
//		}		
	}
	
	public static void main (String args[]) throws URISyntaxException, IOException, GBrowserException {
		
		URL url = new File(System.getProperty("user.home") + "/chipster/Homo_sapiens.GRCh37.66.gtf").toURI().toURL();
		RandomAccessLineDataSource file = new RandomAccessLineDataSource(url, GtfToFeatureConversion.class);
		
		test(file);
		
		url = new URL("http://chipster-filebroker.csc.fi:7060/public/annotations/tmp/Homo_sapiens.GRCh37.66.gtf");
		file = new RandomAccessLineDataSource(url, GtfToFeatureConversion.class);
		
		test(file);

	}

	private static void test(RandomAccessLineDataSource file)
			throws IOException, GBrowserException {
		long t = System.currentTimeMillis();
		
		System.out.println("First 100 lines: ");
		file.setLineReaderPosition(0);
		System.out.println("\t1: " + file.getNextLine());
		
		for (int i = 0; i < 100; i++) {
			file.getNextLine();
		}
		System.out.println("\t100: " + file.getNextLine());
		
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
