package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.net.URL;

import org.apache.commons.io.input.BoundedInputStream;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.UnsortedDataException;
import fi.csc.microarray.util.IOUtils;

public class RandomAccessLineDataSource extends DataSource {

	private Long length = null;
	private BufferedReader lineReader;

	public RandomAccessLineDataSource(URL url, Class<? extends AreaRequestHandler> requestHandler) throws FileNotFoundException, URISyntaxException {
		super(url, requestHandler);
	}

	public RandomAccessLineDataSource(URL urlRoot, String path, Class<? extends AreaRequestHandler> requestHandler)
			throws FileNotFoundException, MalformedURLException, URISyntaxException {
		
		super(urlRoot, path, requestHandler);
	}	
	
	public void setLineReaderPosition(long position) throws IOException, GBrowserException {
		
		if (lineReader != null) {
			IOUtils.closeIfPossible(lineReader);
		}

		if (file != null) {
			FileInputStream in = new FileInputStream(file);
			in.skip(position);
			lineReader = new BufferedReader(new InputStreamReader(in));

		} else {

			HttpURLConnection connection = (HttpURLConnection)url.openConnection();
			InputStream in = connection.getInputStream();
			in.skip(position);
			lineReader = new BufferedReader(new InputStreamReader(in));
			
		}

		if (!lineReader.markSupported()) {
			throw new GBrowserException("File '" + file + "' does't support random access line reading.");
		}
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
		
		while ((line = lineReader.readLine()) != null) {
			lastLine = line;
		}
		
		return lastLine;
	}

	public void close() {
		if (lineReader != null) {
			try {
				lineReader.close();
			} catch (IOException e) {
				//No problem
			}
			lineReader = null;
		}
	}
	
	public long length() throws IOException {
		if (file != null) {
			return file.length();

		} else {
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
		
		System.out.println("First 3 lines: ");
		file.setLineReaderPosition(0);
		System.out.println("\t1: " + file.getNextLine());
		System.out.println("\t2: " + file.getNextLine());
		System.out.println("\t3: " + file.getNextLine());
		
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
