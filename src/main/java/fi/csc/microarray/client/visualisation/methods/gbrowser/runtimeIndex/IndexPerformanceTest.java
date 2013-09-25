package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.util.Collection;
import java.util.Iterator;
import java.util.Random;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;

public class IndexPerformanceTest {
	
	/**
	 * Requires -Xmx2048m
	 * 
	 * 
	 * 
	 * Example output on ThinkPad X220:
	 * 
	 * Memory usage: 0 MB
	 * Init memIndex: 24705 ms
	 * Memory usage: 1528 MB
	 * Init  searchIndex: 10 ms
	 * Memory usage: 1528 MB
	 * Troughput test
	 * memIndex:     	Bytes: 20 MB 	Time: 54 ms 	Bandwidth: 377.79705612747756MB/s
	 * searchIndex:  	Bytes: 20 MB 	Time: 3575 ms 	Bandwidth: 5.706584903743717MB/s
	 * Random access test
	 * memIndex:     Seek 1000, first: 18 ms
	 * Memory usage: 1573 MB
	 * memIndex:    Seek 1000, second: 11 ms
	 * Memory usage: 1576 MB
	 * searchIndex:  Seek 1000, first: 13673 ms
	 * Memory usage: 1587 MB
	 * searchIndex: Seek 1000, second: 13862 ms
	 * Memory usage: 1540 MB
	 * Agreement test
	 * Agreement test done
	 * 
	 * @param args
	 * @throws FileNotFoundException
	 * @throws MalformedURLException
	 * @throws IOException
	 * @throws URISyntaxException
	 * @throws GBrowserException
	 */
	public static void main(String[] args) throws FileNotFoundException, MalformedURLException, IOException, URISyntaxException, GBrowserException {
		
		printMemoryUsage();
		
		long t = System.currentTimeMillis();

		//File has to be sorted
		File file = new File(System.getProperty("user.home") + "/chipster/Homo_sapiens.GRCh37.69-sort.gtf");
		
		InMemoryIndex memIndex;
		BinarySearchIndex searchIndex;
		
		DataUrl dataUrl = new DataUrl(file);

		memIndex = new InMemoryIndex(new LineDataSource(dataUrl), new GtfLineParser());
		
		System.out.println("Init memIndex: " + (System.currentTimeMillis() - t)  + " ms");
		t = System.currentTimeMillis();
		printMemoryUsage();;
		
		searchIndex = new BinarySearchIndex(new RandomAccessLineDataSource(dataUrl), new GtfLineParser());
		
		System.out.println("Init  searchIndex: " + (System.currentTimeMillis() - t)  + " ms");
		t = System.currentTimeMillis();
		printMemoryUsage();		

		//First 100 MB of chr1		
		Region chr1 = new Region(0l, 100*1000*1000l, new Chromosome("1"));
				
		System.out.println("Troughput test");
		troughputTest("memIndex:    ", memIndex, chr1);
		troughputTest("searchIndex: ", searchIndex, chr1);
		
		System.out.println("Random access test");
		seekTest("memIndex:    ", memIndex);
		seekTest("searchIndex: ", searchIndex);
		
		System.out.println("Agreement test");
		agreementTest(memIndex, searchIndex);
		
	}

	private static void troughputTest(String description, Index index, Region chr1) throws IOException, GBrowserException {
		long t = System.currentTimeMillis();

		long bytes = 0;
		
		for (String line : index.getFileLines(chr1).values()) {
			bytes += line.length() + 1;
		}
		
		long dt = (System.currentTimeMillis() - t); 
		
		System.out.println(description + " \tBytes: " + bytes/1024/1024 + " MB \tTime: " + dt + " ms \tBandwidth: " + (bytes/1024.0/1024.0)/(dt/1000.0) + "MB/s");
	}
	
	private static void seekTest(String description, Index index) throws IOException, GBrowserException {
		long t = System.currentTimeMillis();
		
		randomSeeks(index);
		
		System.out.println(description + " Seek 1000, first: " + (System.currentTimeMillis() - t)  + " ms");
		t = System.currentTimeMillis();
		printMemoryUsage();
		
		randomSeeks(index);
		
		System.out.println(description + "Seek 1000, second: " + (System.currentTimeMillis() - t)  + " ms");
		t = System.currentTimeMillis();
		printMemoryUsage();
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
			
			if (!testCompare(lines1, lines2)) {
				System.err.println("\tRequest region: " + region);
			} 
		}		
		System.out.println("Agreement test done");
	}
	
	private static boolean testCompare(Collection<String> lines1, Collection<String> lines2) {
		if (lines1.size() != lines2.size()) {
			System.err.println("Unequal line count: " + lines1.size() + ", " + lines2.size());

		}

		Iterator<String> iter1 = lines1.iterator();
		Iterator<String> iter2 = lines2.iterator();

		while (iter1.hasNext()) {

			String line1 = iter1.next();
			String line2 = iter2.next();

			if (!line1.equals(line2)) {
				System.err.println("Unequal lines");
				System.err.println("\t" + line1);
				System.err.println("\t" + line2);

				return false;
			}
		}

		return true;
	}

	private static void printMemoryUsage() {
		Runtime rt = Runtime.getRuntime();
		System.out.println("Memory usage: " + ((rt.totalMemory() - rt.freeMemory()) / 1024 / 1024) + " MB");
	}
}