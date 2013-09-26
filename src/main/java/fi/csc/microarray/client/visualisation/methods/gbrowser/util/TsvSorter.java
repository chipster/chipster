package fi.csc.microarray.client.visualisation.methods.gbrowser.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;

import org.springframework.util.StringUtils;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.AbstractTsvLineParser;
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.Strings;

public class TsvSorter {

	private int chrCol;
	private int bpCol;
	private ChromosomeNormaliser chromosomeNormaliser = new ChromosomeNormaliser() {

		public String normaliseChromosome(String chromosomeName) {

			// Leave prefix as it is
			
			// Remove postfix, if present
			String SEPARATOR = ".";
			if (chromosomeName.contains(SEPARATOR)) {
				chromosomeName = chromosomeName.substring(0, chromosomeName.indexOf(SEPARATOR));
			}
			
			return chromosomeName;
		}
	};
	
	private AbstractTsvLineParser parser;
	
	public void sort(File in, File out, int chrColumn, int startColumn) throws Exception {
		this.chrCol = chrColumn;
		this.bpCol = startColumn;
		externalSort(in, out);
	}
	
	public void sort(File in, File out, int chrColumn, int startColumn, AbstractTsvLineParser parser) throws Exception {
		this.parser = parser;
		sort(in, out, chrColumn, startColumn);
	}	

	private class Row extends BpCoord {

		public String line;

		public Row(String line) {
			super(null, null);

			this.line = line;
			String[] splitted = line.split("\t");
			String chrStr = splitted.length > chrCol ? splitted[chrCol] : "";
			
			// If chromosome name exists, normalise it
			if (!chrStr.isEmpty()) {
				chrStr = chromosomeNormaliser.normaliseChromosome(chrStr);
				splitted[chrCol] = chrStr;
				this.line = Strings.delimit(Arrays.asList(splitted), "\t"); // replace back to raw line
			}
			String bpStr = splitted.length > bpCol ? splitted[bpCol] : "";

			chr = new Chromosome(chrStr);

			if (bpStr.isEmpty()) {
				bp = -1l;
			} else {
				bp = Long.parseLong(bpStr);
			}
		}
	}

	/**
	 * Based on http://www.codeodor.com/index.cfm/2007/5/14/Re-Sorting-really-BIG-files---the-Java-source-code/1208
	 * 
	 * @param infile
	 * @param outfile
	 * @throws IOException
	 * @throws GBrowserException
	 */
	private void externalSort(File infile, File outfile) throws IOException, GBrowserException {
		
		// Start reading
		BufferedReader initReader = new BufferedReader(new FileReader(infile));
		
		// Read header, if exists
		String header = "";

		if (parser != null) {
			
			String line;
			while ((line = initReader.readLine()) != null) {
				
				parser.setLine(line);
				if (parser.isContentLine()) {
					break;

				} else {
					header += line + "\n";
				}			
			}
			
			//First content line is already consumed. Open file again, but do not read content yet.
			initReader.close();
			initReader = new BufferedReader(new FileReader(infile));
			
			for (int i = 0; i < StringUtils.countOccurrencesOf(header, "\n"); i++) {
				initReader.readLine();
			}					
		}

		// Create and sort chunks
		ArrayList<Row> rowBatch = new ArrayList<Row>(500000);
		boolean quit = false;
		int numFiles = 0;

		while (!quit) {

			// showProgress("Reading...");

			// limit chunks to 100MB
			int size = 0;
			while (size < 100000000) {
				// while (size < 10000000) {
				String line = initReader.readLine();

				if (line == null) {
					quit = true;
					break;
				}

				rowBatch.add(new Row(line));
				size += line.length();
			}

			// showProgress("Sorting...");

			// Use Java's sort.
			Collections.sort(rowBatch);

			// showProgress("Writing...");

			// write to disk
			FileWriter fw = new FileWriter(infile + "_chunk" + numFiles);
			BufferedWriter bw = new BufferedWriter(fw);
			for (int i = 0; i < rowBatch.size(); i++) {
				bw.append(rowBatch.get(i).line + "\n");
			}
			bw.close();
			numFiles++;
			rowBatch.clear();
		}

		// showProgress("Merging...");

		mergeFiles(infile.getAbsolutePath(), outfile, numFiles, header);

		// showProgress("DONE");

		initReader.close();
	}

	private void mergeFiles(String inputFilePath, File outputFilePath, int numChunkFiles, String header) throws IOException, GBrowserException {

		
		// Initialise
		ArrayList<BufferedReader> mergefbr = new ArrayList<BufferedReader>();
		ArrayList<Row> filerows = new ArrayList<Row>();
		FileWriter fw = new FileWriter(outputFilePath);
		BufferedWriter bw = new BufferedWriter(fw);

		try {
			// Write header, if needed
			if (!header.isEmpty()) {
				bw.append(header);	
			}		

			// Merge chunks
			boolean someFileStillHasRows = false;
			for (int i = 0; i < numChunkFiles; i++) {
				mergefbr.add(new BufferedReader(new FileReader(inputFilePath + "_chunk" + i)));

				// get the first row
				String line = mergefbr.get(i).readLine();
				if (line != null) {
					filerows.add(new Row(line));
					someFileStillHasRows = true;
				} else {
					filerows.add(null);
				}
			}

			Row row;
			while (someFileStillHasRows) {
				Row min;
				int minIndex = 0;

				row = filerows.get(0);
				if (row != null) {
					min = row;
					minIndex = 0;
				} else {
					min = null;
					minIndex = -1;
				}

				// check which one is min
				for (int i = 1; i < filerows.size(); i++) {
					row = filerows.get(i);
					if (min != null) {

						if (row != null && (row.compareTo(min) < 0)) {
							minIndex = i;
							min = filerows.get(i);
						}
					} else {
						if (row != null) {
							min = row;
							minIndex = i;
						}
					}
				}

				if (minIndex < 0) {
					someFileStillHasRows = false;
				} else {
					// write to the sorted file
					bw.append(filerows.get(minIndex).line + "\n");

					// get another row from the file that had the min
					String line = mergefbr.get(minIndex).readLine();
					if (line != null) {
						filerows.set(minIndex, new Row(line));
					} else {
						filerows.set(minIndex, null);
					}
				}
				// check if one still has rows
				for (int i = 0; i < filerows.size(); i++) {

					someFileStillHasRows = false;
					if (filerows.get(i) != null) {
						if (minIndex < 0) {
							
							bw.close();
							
							throw new GBrowserException("Error in sorting: " + "mindex lt 0 and found row not null" + filerows.get(i));
						}
						someFileStillHasRows = true;
						break;
					}
				}

				// check the actual files one more time
				if (!someFileStillHasRows) {
					// write the last one not covered above
					for (int i = 0; i < filerows.size(); i++) {
						if (filerows.get(i) == null) {
							String line = mergefbr.get(i).readLine();
							if (line != null) {

								someFileStillHasRows = true;
								filerows.set(i, new Row(line));
							}
						}
					}
				}
			}
		} finally {

			// close all the files			
			IOUtils.closeIfPossible(bw);
			IOUtils.closeIfPossible(fw);
		}

		for (int i = 0; i < mergefbr.size(); i++)
			mergefbr.get(i).close();

		// Delete all of the chunk files.
		for (int i = 0; i < numChunkFiles; i++) {
			File f = new File(inputFilePath + "_chunk" + i);
			f.delete();
		}
	}

	public static void main(String[] args) throws Exception {

		try {
		File in = new File(args[0]);
		File out = new File(args[1]);
		int chr = Integer.parseInt(args[2]);
		int start = Integer.parseInt(args[3]);
		
		new TsvSorter().sort(in, out, chr, start);
		
		} catch (Exception e) {
			e.printStackTrace();
						
			System.out.println(
					"usage: \n" +
					"  TsvSorter <file-in> <file-out> <chr-column> <start-position-column>\n" +
					"  Column indexes start from 0\n\n" +
					"example:\n " +
					"  java -cp chipster-2.7.1.jar fi.csc.microarray.client.visualisation.methods.gbrowser.util.TsvSorter Homo_sapiens.GRCh37.70.gtf Homo_sapiens.GRCh37.70-sort.gtf 0 3");
		}				
	}
}