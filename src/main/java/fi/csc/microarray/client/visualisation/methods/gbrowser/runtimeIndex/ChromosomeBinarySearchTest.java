package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.Iterator;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;

public class ChromosomeBinarySearchTest {
	
	public static void main(String[] args) throws IOException, GBrowserException, URISyntaxException {

		File file = getTestFile();
		
		DataUrl dataUrl = new DataUrl(file);
		
		ChromosomeBinarySearch search = new ChromosomeBinarySearch(dataUrl, new GtfLineParser());
		TreeSet<Chromosome> chrs = search.getChromosomes();
		
		Iterator<Chromosome> chrIter = chrs.iterator();
		
		for (int i = 0; i < 999; i++) {
			Chromosome chr = chrIter.next();
			if (!("" + i).equals(chr.toNormalisedString())) {
				System.err.println("Missing chromosome: " + i);
			}			
		}
		
		file.delete();
	}
	
	
	/**
	 * Create a temp file with 1000 rows:
	 * 
	 * chr0	Source000	Feature000	000	000	Quality000	Strand000	Frame000	Metadata000	
	 * chr1	Source001	Feature001	001	001	Quality001	Strand001	Frame001	Metadata001	
	 * chr2	Source002	Feature002	002	002	Quality002	Strand002	Frame002	Metadata002	
	 * ...
	 * chr999	Source999	Feature999	999	999	Quality999	Strand999	Frame999	Metadata999	
	 * 
	 * @return
	 * @throws IOException
	 */
	public static File getTestFile() throws IOException {
		//Generate test file
		File testFile = File.createTempFile("ChromosomeBinarySearchTest-file", ".txt");
		
		BufferedWriter writer = new BufferedWriter(new FileWriter(testFile));
		
		String[] cols = new String[] { "Source", "Feature", "", "", "Quality", "Strand", "Frame", "Metadata" };
		
		int TEST_FILE_ROWS = 1000;		
		
		for (int i = 0; i < TEST_FILE_ROWS ; i++) {
			String row = String.format("%03d", i);
						
			String line = "";
			
			line += "chr" + i + "\t";
			
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
