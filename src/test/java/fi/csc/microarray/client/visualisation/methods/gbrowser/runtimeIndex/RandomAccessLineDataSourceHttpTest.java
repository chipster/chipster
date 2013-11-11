package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.MalformedURLException;
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

		URL httpUrl = getRemoteGtfUrl();
		URL fileUrl = getLocalGtfUrl();
		
		DataUrl fileDataUrl = new DataUrl(fileUrl, "file");
		RandomAccessLineReader fileReader = new RandomAccessLineReader(fileDataUrl);		
		
		DataUrl httpDataUrl = new DataUrl(httpUrl, "http");
		RandomAccessLineReader httpReader = new RandomAccessLineReader(httpDataUrl);
		
		iterativeTest(fileReader, httpReader);
		
		removeLocalGtfFile(fileUrl);
	}
	
	public static URL getRemoteGtfUrl() throws MalformedURLException {
		/* If this files gets lost by accident, download a new one from ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/
		 * and pick every 1000th line of it to make it smaller:
		 * cat Homo_sapiens.GRCh37.70-sort.gtf | sed -n '0~1000p' > Homo_sapiens.GRCh37.70-500k.gtf
		 */

		return new URL("http://nic.funet.fi/pub/sci/molbio/chipster/devel/junit-test-data/Homo_sapiens.GRCh37.70-500k.gtf");
		//return new URL("http://nic.funet.fi/pub/sci/molbio/chipster/devel/junit-test-data/Homo_sapiens.GRCh37.70-sort.gtf");		
	}
	
	public static URL getLocalGtfUrl() throws MalformedURLException, IOException {
		File file = download(getRemoteGtfUrl());
		return file.toURI().toURL();
	}
	
	public static void removeLocalGtfFile(URL url) throws URISyntaxException {
		new File(url.toURI()).delete();
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
