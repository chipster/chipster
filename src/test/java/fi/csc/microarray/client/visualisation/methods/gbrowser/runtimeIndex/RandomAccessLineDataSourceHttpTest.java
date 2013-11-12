package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.util.Random;

import org.junit.Assert;
import org.junit.Test;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;


public class RandomAccessLineDataSourceHttpTest {
	@Test
	public void run () throws URISyntaxException, IOException, GBrowserException {

		URL httpUrl = getRemoteGtfUrl();
		URL fileUrl = getLocalGtfUrl();
		URL bigHttpUrl = getBigRemoteGtfUrl();
		
		DataUrl fileDataUrl = new DataUrl(fileUrl, "file");
		RandomAccessLineReader fileReader = new RandomAccessLineReader(fileDataUrl);		
		
		DataUrl httpDataUrl = new DataUrl(httpUrl, "http");
		RandomAccessLineReader httpReader = new RandomAccessLineReader(httpDataUrl);
		
		iterativeTest(fileReader, httpReader);
		
		removeLocalGtfFile(fileUrl);
		
		DataUrl bigHttpDataUrl = new DataUrl(bigHttpUrl, "http");
		RandomAccessLineReader bigHttpReader = new RandomAccessLineReader(bigHttpDataUrl);
		testRemoteRandomAccess(bigHttpReader);
	}

	public static URL getRemoteGtfUrl() throws MalformedURLException {
		/* If this files gets lost by accident, download a new one from ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/,
		 * sort it in Chipster and pick every 1000th line of it to make it smaller:
		 * cat Homo_sapiens.GRCh37.70-sort.gtf | sed -n '0~1000p' > Homo_sapiens.GRCh37.70-500k.gtf
		 */

		return new URL("http://nic.funet.fi/pub/sci/molbio/chipster/devel/junit-test-data/Homo_sapiens.GRCh37.70-500k.gtf");		
	}
	
	public static URL getBigRemoteGtfUrl() throws MalformedURLException {
		/* If this files gets lost by accident, download a new one from ftp://ftp.ensembl.org/pub/current_gtf/homo_sapiens/,
		 * sort it in Chipster.
		 */

		return new URL("http://nic.funet.fi/pub/sci/molbio/chipster/devel/junit-test-data/Homo_sapiens.GRCh37.70-sort.gtf");		
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
			
			for (int i = 0; i < 500; i++) {
				String fileLine = file.readLine();
				String httpLine = http.readLine();

				Assert.assertEquals(fileLine, httpLine);
			}
		}				
	}
	
	/**
	 * There should be no difference between random access reading in the beginning of the file
	 * and in the end of the file.
	 * 
	 * @param reader
	 * @throws IOException
	 * @throws GBrowserException
	 */
	private static void testRemoteRandomAccess(RandomAccessLineReader reader) throws IOException, GBrowserException {
						
		Random rand = new Random();
		
		long size = reader.length();	
		final int MB = 1024*1024;
		final long SAMPLES = 10;
		
		long t = System.currentTimeMillis();		
		for (int i = 0; i < SAMPLES; i++) {
			reader.setPosition(rand.nextInt(MB));
			reader.readLine();
		}		
		long headSeek = (System.currentTimeMillis() - t) / SAMPLES;
		
		t = System.currentTimeMillis();
		
		for (int i = 0; i < SAMPLES; i++) {
			reader.setPosition(size - rand.nextInt(MB));
			reader.readLine();
		}		
		long tailSeek = (System.currentTimeMillis() - t) / SAMPLES;					
		
		Assert.assertFalse("Random access read in the end of the file is too slow (" + tailSeek + "ms) " +
				"in comparison to reading in the beginning of the file (" + headSeek + "ms)", 
				tailSeek > headSeek * 2);
	}
}
