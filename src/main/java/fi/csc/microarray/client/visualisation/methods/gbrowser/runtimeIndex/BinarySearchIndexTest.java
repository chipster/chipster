package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.util.Collection;
import java.util.SortedMap;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.IndexKey;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.UnsortedDataException;

public class BinarySearchIndexTest {
		
	public static void main(String[] args) throws FileNotFoundException, MalformedURLException, IOException, URISyntaxException, GBrowserException {
		
		//Test without header
		File testFile = getTestFile(0);		
		runTests(testFile);		
		testFile.delete();
		
		//Test with a 1 line header
		testFile = getTestFile(1);		
		runTests(testFile);		
		testFile.delete();
		
		//Test with a 20 line header
		testFile = getTestFile(20);		
		runTests(testFile);		
		testFile.delete();
	}

	private static void runTests(File file) throws IOException,
			GBrowserException, UnsortedDataException, FileNotFoundException,
			URISyntaxException, MalformedURLException {
		
		DataUrl dataUrl = new DataUrl(file);
		
		Index index = new BinarySearchIndex(new RandomAccessLineDataSource(dataUrl), new GtfLineParser());		
		//Index index = new InMemoryIndex(new LineDataSource(testFile.toURI().toURL(), null), new StackGtfParser());
		
		//Empty region
		Region region = new Region(1l, 1l, new Chromosome("chr1"));
		Collection<String> lines = index.getFileLines(region).values();
		System.out.println(lines.size() == 0);
		
		//Illegal region		
		boolean illegalRegionException = false;
		try {
			region = new Region(2l, 1l, new Chromosome("chr1"));
			index.getFileLines(region);
		} catch (IllegalArgumentException e) {
			illegalRegionException = true;
		}
		System.out.println(illegalRegionException);
				
		//First line
		region = new Region(0l, 1l, new Chromosome("chr1"));
		lines = index.getFileLines(region).values();
		System.out.println(lines.size() == 1);
		System.out.println("chr1\tsource\tfeature\t0\t0\tQuality0000\tStrand0000\tFrame0000\tMetadata0000\t".equals(lines.iterator().next()));
		
		//Second line
		region = new Region(1l, 2l, new Chromosome("chr1"));
		lines = index.getFileLines(region).values();
		System.out.println(lines.size() == 1);
		System.out.println("chr1\tsource\tfeature\t1\t1\tQuality0001\tStrand0001\tFrame0001\tMetadata0001\t".equals(lines.iterator().next()));
		
		//Several lines
		region = new Region(1000l, 2000l, new Chromosome("chr1"));
		lines = index.getFileLines(region).values();
		System.out.println(lines.size() == 1000);
			
		//Request equals index entry (there are 2000 lines with start position 4000
		//Index entry is somewhere around line 5000, but this should still give also 1000 equal lines (in sort order) before it
		region = new Region(4000l, 6000l, new Chromosome("chr1"));
		lines = index.getFileLines(region).values();
		System.out.println(lines.size() == 2000);
				
		//Last but one
		region = new Region(9998l, 9999l, new Chromosome("chr1"));
		lines = index.getFileLines(region).values();
		System.out.println(lines.size() == 1);
		System.out.println("chr1\tsource\tfeature\t9998\t9998\tQuality9998\tStrand9998\tFrame9998\tMetadata9998\t".equals(lines.iterator().next()));
				
		//Last line
		region = new Region(9999l, 10000l, new Chromosome("chr1"));
		lines = index.getFileLines(region).values();
		System.out.println(lines.size() == 1);
		System.out.println("chr1\tsource\tfeature\t9999\t9999\tQuality9999\tStrand9999\tFrame9999\tMetadata9999\t".equals(lines.iterator().next()));
		
		//Region greater than last line
		region = new Region(10000l, 20000l, new Chromosome("chr1"));
		lines = index.getFileLines(region).values();
		System.out.println(lines.size() == 0);
		
		//Test line identifiers on the first line of the file, last line and in between.		
		testLineIdentifiers(index, 0);
		testLineIdentifiers(index, 1000);
		testLineIdentifiers(index, 9999);
	}

	private static void testLineIdentifiers(Index index, long startPosition) throws IOException,
			GBrowserException {
		
		//Line identifier must be same regardless of the line position in the request region 
		
		//First line of request region
		Region region;
		region = new Region(startPosition, startPosition + 100, new Chromosome("chr1"));
		TreeMap<IndexKey, String> lineMap = index.getFileLines(region);
		IndexKey id1 = lineMap.firstKey();
		
		//Last line  of request region
		region = new Region(startPosition - 100, startPosition + 1, new Chromosome("chr1"));
		lineMap = index.getFileLines(region);
		IndexKey id2 = lineMap.lastKey();
		
		//In the middle of request region
		region = new Region(startPosition - 100, startPosition + 100, new Chromosome("chr1"));
		lineMap = index.getFileLines(region);
		IndexKey fromKey = new IndexKey(new BpCoord(startPosition, new Chromosome("chr1")), 0);
		IndexKey toKey = new IndexKey(new BpCoord(startPosition, new Chromosome("chr1")), Long.MAX_VALUE);
		SortedMap<IndexKey, String> subMap = lineMap.subMap(fromKey, toKey);		
		IndexKey id3 = subMap.firstKey();
		
		System.out.println(id1.equals(id2) && id2.equals(id3));
	}

//	private static void print(List<String> lines) {
//		System.out.println(lines.size() + " line(s):");
//		if (lines.size() >= 1) {
//			System.out.println(lines.get(0));
//		}
//		if (lines.size() >= 2) {
//			System.out.println(lines.get(1));
//		}
//		if (lines.size() >= 3) {
//			System.out.println("...");
//			System.out.println(lines.get(lines.size() - 1));
//		}
//	}
	
	/**
	 * Create a temp file with 10000 rows:
	 * 
	 * chr1	source	feature	0	0	Quality0000	Strand0000	Frame0000	Metadata0000	
	 * chr1	source	feature	1	1	Quality0001	Strand0001	Frame0001	Metadata0001	
	 * chr1	source	feature	2	2	Quality0002	Strand0002	Frame0002	Metadata0002	
	 * ...
	 * chr1	source	feature	9999	9999	Quality9999	Strand9999	Frame9999	Metadata9999
	 * @param headerLineCount 
	 * 
	 * @return
	 * @throws IOException
	 */
	public static File getTestFile(int headerLineCount) throws IOException {
		//Generate test file
		File testFile = File.createTempFile("RandomAccessBinarySearchTest-file", ".txt");		
		BufferedWriter writer = new BufferedWriter(new FileWriter(testFile));
		
		//Header
		for (int i = 0; i < headerLineCount; i++) {
			writer.write("#header line " + i);
			writer.newLine();
		}
		
		String[] cols = new String[] { "Quality", "Strand", "Frame", "Metadata" };
		
		int TEST_FILE_ROWS = 10000;
		
		int startPosition = 0;
		
		for (int i = 0; i < TEST_FILE_ROWS ; i++) {
			String row = String.format("%04d", i);
						
			String line = "";
			
			//Keep start position fixed around first index entry
			if (i <= 4000 || i >= 6000) {
				startPosition = i;
			}
			
			line += "chr1\tsource\tfeature\t";
			line += startPosition + "\t";
			line += startPosition + "\t";
			
			for (String col : cols) {
				line += col + row + "\t";				
			}
			
			//System.out.println(line);
			writer.write(line);
			writer.newLine();
		}
		writer.flush();
		writer.close();
		return testFile;
	}
}