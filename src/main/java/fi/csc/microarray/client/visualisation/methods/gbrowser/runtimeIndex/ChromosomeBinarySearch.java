package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.UnsortedDataException;

/**
 * Sorted files make it possible to find quickly requested region, but there is no way 
 * to know what to request without information about chromosome names. This class searches 
 * chromosme names from sorted file relatively quickly with binary search.
 * 
 * @author klemela
 */
public class ChromosomeBinarySearch {

	//Read through smaller regions than this
	private static final long ITERATIVE_SEARCH_LIMIT = 8*1024;
	
	private RandomAccessLineDataSource file;
	
	//Set of chromosome names, prevents duplicates
	private TreeSet<Chromosome> chrSet = new TreeSet<Chromosome>();

	private LineParser parser;

	public ChromosomeBinarySearch(DataUrl data, LineParser parser)
			throws IOException, GBrowserException, URISyntaxException, UnsortedDataException {
		
		this.file = new RandomAccessLineDataSource(data);					
		this.parser = parser;
	}
	
	public TreeSet<Chromosome> getChromosomes() throws IOException, GBrowserException, UnsortedDataException {
		
		//Only for sorting check
		BinarySearchIndex index = new BinarySearchIndex(file, parser);
		index.checkSorting();
		
		//Do the work
		searchChromosomeChange(0, Long.MAX_VALUE, null, null);
		
		return chrSet;
	}

	private Chromosome getChr(long pos) throws IOException, GBrowserException {
				
		if (pos != 0) {
			file.setLineReaderPosition(pos);			
			//Usually partial line
			file.getNextLine();
		} else {
			file.setLineReaderPosition(0);
		}
		
		String line = file.getNextLine();
		
		if (line != null) {
			parser.setLine(line);
			if (parser.isContentLine()) {
				return parser.getRegion().start.chr;
			} else {
				//Header or comment line, read the next one
				return getChr(pos + line.length() + 1); 
			}
		} else {
			return null;
		}
	}
	
	/**
	 * Search chromosome change recursively between the given file positions. Corresponding
	 * chromosomes can be supplied to avoid unnecessary file reading, or null if they are not known yet.
	 * 
	 * @param pos1
	 * @param pos2
	 * @param chr1
	 * @param chr2
	 * @throws IOException
	 * @throws GBrowserException
	 */
	private void searchChromosomeChange(long pos1, long pos2, Chromosome chr1, Chromosome chr2) throws IOException, GBrowserException {
		
		if (pos1 >= pos2) {
			return;
		}
		
		if (chr1 == null) {
			chr1 = getChr(pos1);
		}
		
		//Special case in search start
		if (pos2 == Long.MAX_VALUE) {
			
			String lastLine = file.getLastLine();
			parser.setLine(lastLine);
			chr2 = parser.getRegion().start.chr;
			//Minus one to point to preceding new line character
			pos2 = file.length() - lastLine.length() - 1;
			
		} else {
			
			if (chr2 == null) {
				chr2 = getChr(pos2);
			}
		}
				
		chrSet.add(chr1);		
		chrSet.add(chr2);
		
		if (chr1.equals(chr2)) {
			//There is no change
		} else {
			
			if (pos2 - pos1 < ITERATIVE_SEARCH_LIMIT) {

				linearSearch(pos1, pos2, chr1, chr2);
				
			} else {
				
				long centerPos = (pos1 + pos2) / 2;		
				
				//This is same for both calls, so its more efficient to read it here
				Chromosome centerChr = getChr(centerPos);
				
				searchChromosomeChange(pos1, centerPos, chr1, centerChr);
				searchChromosomeChange(centerPos, pos2, centerChr, chr2);
			}						
		}
	}

	/**
	 * Search chromosome changes line-by-line, until a chr2 is found.
	 * 
	 * @param pos1
	 * @param pos2
	 * @param chr1
	 * @param chr2
	 * @throws IOException
	 * @throws GBrowserException
	 */
	private void linearSearch(long pos1, long pos2, Chromosome chr1,
			Chromosome chr2) throws IOException, GBrowserException {
			
		long pos = pos1;

		file.setLineReaderPosition(pos);

		//Usually partial line
		if (pos != 0) {
			pos += file.getNextLine().length() + 1; //new line character
		}

		while (pos < pos2) {
			String line = file.getNextLine();
			parser.setLine(line);

			if (!parser.isContentLine()) {
				//header or comment
				continue;
			}
			Chromosome searchChr = parser.getRegion().start.chr;
			pos += line.length() + 1; //new line character

			if (!chr1.equals(searchChr)) {

				//Found new chr 
				chrSet.add(searchChr);
				
				if (chr2.equals(searchChr)) {

					//No more chromosome changes
					break;					
				}								
			}
		}								
	}
	
	/**
	 * Example results for 500 MB gtf file with 200 chromosomes:
	 * 
	 * File: 295 ms
	 * Http: 2386 ms (ping <1ms)
	 * 
	 * @param args
	 * @throws IOException
	 * @throws GBrowserException
	 * @throws URISyntaxException
	 */
	public static void main(String[] args) throws IOException, GBrowserException, URISyntaxException {
		String fileString = System.getProperty("user.home") + "/chipster/Homo_sapiens.GRCh37.69-sort.gtf";
		
		URL fileUrl = new File(fileString).toURI().toURL();			
		printTestTime(fileUrl, "File: ");
		
		URL httpUrl = new URL("http://chipster-filebroker.csc.fi:7060/public/annotations/tmp/Homo_sapiens.GRCh37.69-sort.gtf");
		
		printTestTime(httpUrl, "Http: ");
	}

	private static void printTestTime(URL url, String description) throws IOException,
			GBrowserException, URISyntaxException {
		ChromosomeBinarySearch index = new ChromosomeBinarySearch(new DataUrl(url, url.getPath()), new GtfLineParser());

		long t = System.currentTimeMillis();
		
		@SuppressWarnings("unused")
		TreeSet<Chromosome> chrs = index.getChromosomes();
		
		System.out.println(description + (System.currentTimeMillis() - t) + " ms");
		
//		for (Chromosome chr : chrs) {
//			System.out.println(chr);
//		}
	}
}
