package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;

import org.junit.Assert;
import org.junit.Test;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;


public class RandomAccessLineDataSourceHttpTest {
	@Test
	public void run () throws URISyntaxException, IOException, GBrowserException {
				
		URL httpUrl = new URL("http://nic.funet.fi/pub/sci/molbio/chipster/devel/junit-test-data/Homo_sapiens.GRCh37.70-500k.gtf");
		//URL httpUrl = new URL("http://nic.funet.fi/pub/sci/molbio/chipster/devel/junit-test-data/Homo_sapiens.GRCh37.70-sort.gtf");
		
		File file = download(httpUrl);
		URL fileUrl = file.toURI().toURL();
		
		DataUrl fileDataUrl = new DataUrl(fileUrl, "file");
		RandomAccessLineReader fileReader = new RandomAccessLineReader(fileDataUrl);		
		
		DataUrl httpDataUrl = new DataUrl(httpUrl, "http");
		RandomAccessLineReader httpReader = new RandomAccessLineReader(httpDataUrl);
		
		iterativeTest(fileReader, httpReader);
		
		file.delete();
	}

	private static File download(URL httpUrl) throws IOException {
		ReadableByteChannel rbc = Channels.newChannel(httpUrl.openStream());
		File out = File.createTempFile("RandomAccessLineDataSourceHttpTest", "");
		FileOutputStream fos = new FileOutputStream(out);
		fos.getChannel().transferFrom(rbc, 0, Long.MAX_VALUE);
		fos.close();
		return out;
	}

	private static void iterativeTest(RandomAccessLineReader file,
			RandomAccessLineReader http) throws IOException, GBrowserException {
		
		for (int j = 0; j < 3; j++){

			file.setPosition(j*100);
			http.setPosition(j*100);
			
			for (int i = 0; i < 1000; i++) {
				String fileLine = file.readLine();
				String httpLine = http.readLine();

				Assert.assertEquals(fileLine, httpLine);
			}
		}				
	}
}
