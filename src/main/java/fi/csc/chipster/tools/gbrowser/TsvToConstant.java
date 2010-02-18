package fi.csc.chipster.tools.gbrowser;



import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

import fi.csc.microarray.exception.MicroarrayException;


public class TsvToConstant {

	public static File convert(File in, File out, int[] fieldLengths) throws MicroarrayException {
		
		//File out = new File(in.getName().replace(".tsv", ".fsf"));
		
		//System.out.print("Converting to constant column length...");
		//Long startTime = System.currentTimeMillis();

		BufferedWriter writer;
		BufferedReader reader;
		try {
			writer = new BufferedWriter(new FileWriter(out));
			reader = new BufferedReader(new FileReader(in));

			int i = 0;
			String line;
			while((line = reader.readLine()) != null){

				String[] cols = line.split("\t");

				StringBuilder constLine = new StringBuilder();

				for (int j = 0; j < fieldLengths.length; j++){
					//System.out.println(j);
					if(cols.length <= j) {
						constLine.append(fillWithSpaces("", fieldLengths[j]));
					} else {
						try {
							constLine.append(fillWithSpaces(cols[j], fieldLengths[j]));
						} catch (StringIndexOutOfBoundsException e) {
							
							throw new MicroarrayException("Error in converting to constant column length file: " +
									"Column in input file is longer than defined column length. \n" +
									"Input file: " + in.getName() + ", on row: " + i + ", column: " + j + "\n" + 
									"Column content: " + cols[j] + ", defined length: " + fieldLengths[j]);									
						}
					}
				}
				constLine.append("\n");

				writer.write(constLine.toString());

//				if( i % 1000000 == 0){
//					System.out.println("" + (int)((line.length() * i) / (float)in.length() * 100) + " % ");
//				}												

				i++;
			}

			writer.flush();
			writer.close();
			reader.close();
			//System.out.println("DONE in " + (System.currentTimeMillis() - startTime) / 1000 + " seconds");
			
			return out;
		} catch (IOException e) {
			e.printStackTrace();
		}
		return null;
	}


	private static String fillWithSpaces(String orig, int length){
		String spaces = "                                                                                 ";
		return orig + spaces.substring(0, length - orig.length());
	}
}
