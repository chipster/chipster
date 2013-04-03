package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;

public class RandomAccessLineDataSourceTest {
	public static void main (String args[]) throws URISyntaxException, IOException, GBrowserException {
		
		//Functionality
		testFunctionality();
		
		//Performance
		URL fileUrl = new File(System.getProperty("user.home") + "/chipster/Homo_sapiens.GRCh37.69-sort.gtf").toURI().toURL();
		RandomAccessLineDataSource file = new RandomAccessLineDataSource(fileUrl);		
		
		URL httpUrl = new URL("http://chipster-filebroker.csc.fi:7060/public/annotations/tmp/Homo_sapiens.GRCh37.69-sort.gtf");
		RandomAccessLineDataSource http = new RandomAccessLineDataSource(httpUrl);
		
		manualTest(file);
		manualTest(http);
		
		iterativeTest(file, http);
	}

	private static void testFunctionality() throws URISyntaxException, IOException, GBrowserException {
				
		File testFile = RandomAccessLineReaderTest.getTestFile();
		List<String> lines = RandomAccessLineReaderTest.getTestReferenceList(testFile);
		
		RandomAccessLineDataSource dataSource = new RandomAccessLineDataSource(testFile.toURI().toURL());
		
		System.out.println(lines.get(lines.size() - 1).equals(dataSource.getLastLine()));
		
		testFile.delete();
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
		file.setLineReaderPosition(100*1024*1024);
		System.out.println("\t1: " + file.getNextLine());
		

		System.out.println(System.currentTimeMillis() - t + " ms ");
		t = System.currentTimeMillis();
		
		System.out.println("Seek to 300MB: ");
		file.setLineReaderPosition(300*1024*1024);
		System.out.println("\t1: " + file.getNextLine());

		System.out.println(System.currentTimeMillis() - t + " ms ");
		t = System.currentTimeMillis();
		
		System.out.println("Last line: ");
		System.out.println("\t1: " + file.getLastLine());
		
		System.out.println(System.currentTimeMillis() - t + " ms ");
		t = System.currentTimeMillis();
	}
}
