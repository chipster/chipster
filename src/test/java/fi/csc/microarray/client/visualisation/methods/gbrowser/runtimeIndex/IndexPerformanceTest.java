package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Random;

import org.junit.Assert;
import org.junit.Test;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;

public class IndexPerformanceTest {
	
	/**
	 * Test indexes' performance. The limits are really generous to tolerate differences in hardware, OS and JVM:
	 * about 10x for the time and about 2x for the memory usage in comparison to Thinkpad X220 with SSD disk, 
	 * Red hat 6 and Oracle Java 1.7.  
	 * 
	 * @param args
	 * @throws FileNotFoundException
	 * @throws MalformedURLException
	 * @throws IOException
	 * @throws URISyntaxException
	 * @throws GBrowserException
	 */
	
	private static long time;
	private static long memoryUsage;
	
	@Test
	public void run() throws FileNotFoundException, MalformedURLException, IOException, URISyntaxException, GBrowserException {		
		
		//init measuring equipment
		checkLimits(null, null);
		
		File file = new File(RandomAccessLineDataSourceHttpTest.getLocalGtfUrl().toURI());
				
		InMemoryIndex memIndex;
		BinarySearchIndex searchIndex;
		
		DataUrl dataUrl = new DataUrl(file);

		memIndex = new InMemoryIndex(new LineDataSource(dataUrl), new GtfLineParser());	
		checkLimits(3000, 12); //X220: 300ms, 6Mb
				
		searchIndex = new BinarySearchIndex(new RandomAccessLineDataSource(dataUrl), new GtfLineParser());
		checkLimits(100, 16); //X220: 9ms, 8Mb
				
		Region chr1 = new Region(0l, Long.MAX_VALUE, new Chromosome("1"));
				
		//Throughput 
		
		/* This is severely limited by the small size of the test file. 
		 * memIndex throughput is here maybe ~5 Mb/s,
		 * while with larger files it can reach for example 400 Mb/s 
		 */
		troughputTest("memIndex:    ", memIndex, chr1);
		checkLimits(3000, 12); //X220: 300ms, 6Mb
				
		troughputTest("searchIndex: ", searchIndex, chr1);
		checkLimits(100, 10); //X220: 11ms, 1Mb
		
		//The second run is faster only because of the JVM optimizations
		randomSeeks(memIndex);
		checkLimits(500, 10); //X220: 55ms, 0Mb
		randomSeeks(memIndex); 
		checkLimits(250, 10); //X220: 22ms, 2Mb
		
		//The second run really should be faster, because the index is largely build already
		randomSeeks(searchIndex);
		checkLimits(10000, 600); //X220: 1100ms, 265Mb
		randomSeeks(searchIndex); 
		checkLimits(5000, 10); //X220: 450ms, 0Mb
		
		agreementTest(memIndex, searchIndex);
		
	}
	
	private static void checkLimits(Integer i, Integer j) {

		long dTime = System.currentTimeMillis() - time;
		long dMemory = getMemoryUsage() - memoryUsage;
		
		//System.out.println(dTime + "ms, " + dMemory + "Mb");
		
		if (i != null) {
			Assert.assertTrue(dTime < i);
		}
		if (j != null) {
			Assert.assertTrue(dMemory < j);
		}
		
		memoryUsage = getMemoryUsage();		
		time = System.currentTimeMillis();
	}

	private static void troughputTest(String description, Index index, Region chr1) throws IOException, GBrowserException {
		long t = System.currentTimeMillis();
		long bytes = 0;
		
		for (String line : index.getFileLines(chr1).values()) {
			bytes += line.length() + 1;
		}
		
		long dt = (System.currentTimeMillis() - t); 
		
		//System.out.println(description + " \tBytes: " + bytes/1024/1024 + " MB \tTime: " + dt + " ms \tBandwidth: " + (bytes/1024.0/1024.0)/(dt/1000.0) + "MB/s");
	}

	private static Random randomSeeks(Index index) throws IOException,
			GBrowserException {
		
		Random rand = new Random();
				
		for (int i = 0; i < 1000; i++) {
			
			Chromosome chr = new Chromosome("" + rand.nextInt(22));
			long start = rand.nextInt(200*1000*1000);
			long end = start + 100*1000;			
			
			Region region = new Region(start, end, chr);
			
			index.getFileLines(region);
		}
		return rand;
	}
	

	private static void agreementTest(InMemoryIndex index1,
			BinarySearchIndex index2) throws IOException, GBrowserException {
		
		Random rand = new Random();

		for (int i = 0; i < 1000; i++) {

			Chromosome chr = new Chromosome("" + rand.nextInt(22));
			long start = rand.nextInt(200*1000*1000);
			long end = start + 100*1000;			

			Region region = new Region(start, end, chr);

			Collection<String> lines1 = index1.getFileLines(region).values();
			Collection<String> lines2 = index2.getFileLines(region).values();

			testCompare(lines1, lines2, "\tRequest region: " + region);
		} 
	}		
	
	private static boolean testCompare(Collection<String> lines1, Collection<String> lines2, String details) {
		if (lines1.size() != lines2.size()) {
			Assert.fail("Unequal line count: " + lines1.size() + ", " + lines2.size() + details);
		}

		Iterator<String> iter1 = lines1.iterator();
		Iterator<String> iter2 = lines2.iterator();

		while (iter1.hasNext()) {

			String line1 = iter1.next();
			String line2 = iter2.next();

			if (!line1.equals(line2)) {
				Assert.fail("Unequal lines" + details + "\n" + line1 + "\n" + line2);
			}
		}

		return true;
	}

	/**
	 * @return memory usage in Mb
	 */
	private static long getMemoryUsage() {
		Runtime rt = Runtime.getRuntime();
		return (rt.totalMemory() - rt.freeMemory()) / 1024 / 1024;
	}
}