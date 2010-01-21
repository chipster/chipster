package fi.csc.microarray.client.visualisation.methods.genomeBrowser.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Collections;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.ElandParser;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.Chromosome;

public class TsvSorter2 {
	
	private static int chrCol;
	private static int bpCol;
	
	public static void main(String[] args) throws Exception {
		
		//Eland
		chrCol = new ElandParser().getFileDefinition().indexOf(ColumnType.CHROMOSOME);
		bpCol = new ElandParser().getFileDefinition().indexOf(ColumnType.BP_START);
		
		//Casava
//		chrCol = 10;
//		bpCol = 12
		
		long s = System.currentTimeMillis();

		externalSort("eland_result.txt", "sorted.txt");

		System.out.println(System.currentTimeMillis() - s);
	}
	
	private static class Row extends BpCoord{
		
		public String line;
		
		public Row(String line){
			super(null, null);
			
			this.line = line;
			String[] splitted = line.split("\t");
			String chrStr = splitted.length > chrCol ? splitted[chrCol] : "";
			String bpStr = splitted.length > bpCol ? splitted[bpCol] : "";
			
			chr = new Chromosome(chrStr.replace("chr", "").replace(".fa", ""));			

			if(bpStr.equals("")){
				bp = -1l;
			} else {
				bp = Long.parseLong(bpStr);
			}
		}
	}
	
	private static long startTime = -1;
	
	private static void showProgress(String message){
		if(startTime != -1){
			System.out.println((System.currentTimeMillis() - startTime) / 1000 + " s");
		}
		startTime = System.currentTimeMillis();
		
		System.out.println(message);
	}
	
	private static void externalSort(String infile, String outfile) {
		try {
			BufferedReader initReader = new BufferedReader(new FileReader(infile));
			ArrayList<Row> rowBatch = new ArrayList<Row>(500000);

			boolean quit = false;
			int numFiles = 0;

			while (!quit) {
				
				showProgress("Reading...");
				
				// limit chunks to 200MB
				int size = 0;
				while (size < 200000000) {
				//while (size < 10000000) {
					String line = initReader.readLine();

					if (line == null) {
						quit = true;
						break;
					}

					rowBatch.add(new Row(line));
					size += line.length();
				}

				showProgress("Sorting...");
				
				// Use Java's sort.
				Collections.sort(rowBatch);

				
				showProgress("Writing...");
				
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
			
			showProgress("Merging...");

			mergeFiles(infile, outfile, numFiles);
			
			showProgress("DONE");

			initReader.close();
		} catch (Exception ex) {
			ex.printStackTrace();
			System.exit(-1);
		}

	}

	private static void mergeFiles(String inputFilePath, String outputFilePath, int numChunkFiles) {
		try {
			ArrayList<BufferedReader> mergefbr = new ArrayList<BufferedReader>();
			ArrayList<Row> filerows = new ArrayList<Row>();
			FileWriter fw = new FileWriter(outputFilePath);
			BufferedWriter bw = new BufferedWriter(fw);

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
							System.out.println("mindex lt 0 and found row not null"
									+ filerows.get(i));
							System.exit(-1);
						}
						someFileStillHasRows = true;
						break;
					}
				}

				// check the actual files one more time
				if (!someFileStillHasRows) {
					//write the last one not covered above
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

			// close all the files
			bw.close();
			fw.close();
			for (int i = 0; i < mergefbr.size(); i++)
				mergefbr.get(i).close();

			// Delete all of the chunk files.
			for (int i = 0; i < numChunkFiles; i++) {
				File f = new File(inputFilePath + "_chunk" + i);
				f.delete();
			}

		} catch (Exception ex) {
			ex.printStackTrace();
			System.exit(-1);
		}
	}
}