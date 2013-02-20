package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;

public class ChromosomeBinarySearch {

	private static final long ITERATIVE_SEARCH_LIMIT = 10*1024;
	
	private RandomAccessLineDataSource file;
	
	private TreeSet<Chromosome> chrSet = new TreeSet<Chromosome>();

	private Parser parser;

	public ChromosomeBinarySearch(URL url, Parser parser)
			throws IOException, GBrowserException, URISyntaxException {
		
		file = new RandomAccessLineDataSource(url, null);	
				
		this.parser = parser;
	}
	
	public TreeSet<Chromosome> getChromosomes() throws IOException, GBrowserException {
		
		BinarySearchIndex index = new BinarySearchIndex(file, parser);
		index.checkSorting();
		
		searchChromosomeChange(0, Long.MAX_VALUE, null, null);
		
		return chrSet;
	}

	private Chromosome getChr(long pos) throws IOException, GBrowserException {
		
		
		//Usually partial line
		if (pos != 0) {
			//Make sure we get the whole line
			file.setLineReaderPosition(pos - 1);
			file.getNextLine();
		} else {
			file.setLineReaderPosition(0);
		}
		
		String line = file.getNextLine();
		
		if (line != null) {
			parser.setLine(line);
			return parser.getRegion().start.chr;
		} else {
			return null;
		}
		
	}
	
	private void searchChromosomeChange(long pos1, long pos2, Chromosome chr1, Chromosome chr2) throws IOException, GBrowserException {
		
		if (chr1 == null) {
			chr1 = getChr(pos1);
		}
		
				
		if (pos2 == Long.MAX_VALUE) {
			
			String lastLine = file.getLastLine();
			parser.setLine(lastLine);
			chr2 = parser.getRegion().start.chr;
			
			pos2 = file.length() - lastLine.length() - 1;
			
		} else {
			
			if (chr2 == null) {
				chr2 = getChr(pos2);
			}
		}
				
		chrSet.add(chr1);
		
		//Null in end of file
		if (chr2 != null) {
			chrSet.add(chr2);
		}
		
		if (chr1.equals(chr2)) {
			//There is no change
		} else {
			
			if (pos2 - pos1 < ITERATIVE_SEARCH_LIMIT) {

				iterativeSearch(pos1, pos2, chr1, chr2);
				
			} else {
				
				long centerPos = (pos1 + pos2) / 2;
				
				searchChromosomeChange(pos1, centerPos, chr1, null);
				searchChromosomeChange(centerPos, pos2, null, chr2);
			}						
		}
	}

	private void iterativeSearch(long pos1, long pos2, Chromosome chr1,
			Chromosome chr2) throws IOException, GBrowserException {
			
		long pos = pos1;

		file.setLineReaderPosition(pos);

		//Usually partial line
		if (pos != 0) {
			pos += file.getNextLine().length() + 1;
		}

		while (pos < pos2) {
			String line = file.getNextLine();
			parser.setLine(line);
			Chromosome searchChr = parser.getRegion().start.chr;
			pos += line.length() + 1;

			if (!chr1.equals(searchChr)) {

				//Found new chr 
				chrSet.add(searchChr);
				

				if (chr2.equals(searchChr)) {

					//There is another change
					searchChromosomeChange(pos, pos2, searchChr, chr2);
				}								
			}
		}								
	}
	
	public static void main(String[] args) throws IOException, GBrowserException, URISyntaxException {
		String fileString = System.getProperty("user.home") + "/chipster/Homo_sapiens.GRCh37.66-sort.gtf";
		//String fileString = System.getProperty("user.home") + "/chipster/cufflinks-gtf/merged-sort.gtf";
		//String fileString = System.getProperty("user.home") + "/chipster/cufflinks-gtf/transcripts-sort.gtf";
		
		URL url = new File(fileString).toURI().toURL();
		
		//URL url = new URL("http://chipster-filebroker.csc.fi:7060/public/annotations/tmp/Homo_sapiens.GRCh37.66.gtf");
			
		ChromosomeBinarySearch index = new ChromosomeBinarySearch(url, new StackGtfParser());

		long t = System.currentTimeMillis();
		
		TreeSet<Chromosome> chrs = index.getChromosomes();
		
		System.out.println((System.currentTimeMillis() - t) + " ms");
		
		for (Chromosome chr : chrs) {
			System.out.println(chr);
		}
	}
}
