package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;


public class RandomAccessLineDataSourceHttpTest {
	public static void main (String args[]) throws URISyntaxException, IOException, GBrowserException {
		
		URL fileUrl = new File(System.getProperty("user.home") + "/chipster/Homo_sapiens.GRCh37.69-sort.gtf").toURI().toURL();
		DataUrl fileDataUrl = new DataUrl(fileUrl, "file");
		RandomAccessLineReader file = new RandomAccessLineReader(fileDataUrl);		
		
		URL httpUrl = new URL("http://chipster-filebroker.csc.fi:7060/public/annotations/tmp/Homo_sapiens.GRCh37.69-sort.gtf");
		DataUrl httpDataUrl = new DataUrl(httpUrl, "http");
		RandomAccessLineReader http = new RandomAccessLineReader(httpDataUrl);
		
		iterativeTest(file, http);
	}

	private static void iterativeTest(RandomAccessLineReader file,
			RandomAccessLineReader http) throws IOException, GBrowserException {
		
		for (int j = 0; j < 3; j++){
			System.out.println(j);
			file.setPosition(j*100);
			http.setPosition(j*100);
			
			for (int i = 0; i < 1000; i++) {
				String fileLine = file.readLine();
				String httpLine = http.readLine();

				testCompare(fileLine, httpLine, "" + i);
			}
		}
		
		System.out.println("Test done");		
	}

	private static void testCompare(String line1, String line2, String description) {
		if (!line1.equals(line2)) {
			System.out.println("Lines are different!" + description);
			System.out.println("\t" + line1);
			System.out.println("\t" + line2);
		}
	}
}
