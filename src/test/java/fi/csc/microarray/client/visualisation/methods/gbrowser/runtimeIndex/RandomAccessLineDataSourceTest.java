package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.List;

import org.junit.Assert;
import org.junit.Test;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;

public class RandomAccessLineDataSourceTest {
	
	@Test
	public void run () throws URISyntaxException, IOException, GBrowserException {
		
		//Functionality
		testLastLine();
				
		//Performance
		URL fileUrl = RandomAccessLineDataSourceHttpTest.getLocalGtfUrl();
		DataUrl fileDataUrl = new DataUrl(fileUrl, "file");
		
		RandomAccessLineDataSource file = new RandomAccessLineDataSource(fileDataUrl);		
		
		URL httpUrl = RandomAccessLineDataSourceHttpTest.getRemoteGtfUrl();
		DataUrl httpDataUrl = new DataUrl(httpUrl, "http");
		RandomAccessLineDataSource http = new RandomAccessLineDataSource(httpDataUrl);
		
		remoteVsLocal(file, http);
		
		RandomAccessLineDataSourceHttpTest.removeLocalGtfFile(fileUrl);
	}

	private static void testLastLine() throws URISyntaxException, IOException, GBrowserException {
				
		File testFile = RandomAccessLineReaderTest.getTestFile();
		DataUrl dataUrl = new DataUrl(testFile);
		List<String> lines = RandomAccessLineReaderTest.getTestReferenceList(testFile);
		
		RandomAccessLineDataSource dataSource = new RandomAccessLineDataSource(dataUrl);
		
		Assert.assertEquals(lines.get(lines.size() - 1), dataSource.getLastLine());
		
		testFile.delete();
	}

	private static void remoteVsLocal(RandomAccessLineDataSource file,
			RandomAccessLineDataSource http) throws IOException, GBrowserException {
		
		testCompare(file.getLastLine(), http.getLastLine(), "lastLine");		
	}

	private static void testCompare(String line1, String line2, String description) {
		
		if (!line1.equals(line2)) {
			Assert.fail("Lines are different!" + description);
			System.out.println("\t" + line1);
			System.out.println("\t" + line2);
		}
	}
}
