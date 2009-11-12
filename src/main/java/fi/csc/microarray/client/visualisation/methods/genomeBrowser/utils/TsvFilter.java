package fi.csc.microarray.client.visualisation.methods.genomeBrowser.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

public class TsvFilter {
	public static void main(String[] args) {

		System.out.println("File filtering started");
		Long startTime = System.currentTimeMillis();

		//File in = new File("eland_result.txt");
		File in = new File("chr1_filtered.txt");
		File out = new File("u_filtered.txt");		
		int col = 2;
		
		List<String> whatToFind = Arrays.asList(new String[]{
			"U0", "U1", "U2"
		});

		BufferedWriter writer;
		BufferedReader reader;
		try {
			writer = new BufferedWriter(new FileWriter(out));
			reader = new BufferedReader(new FileReader(in));

			int i = 0;
			int accepted = 0;
			String line;
			while((line = reader.readLine()) != null){
				
				//System.out.println(line);

				String[] cols = line.split("\t");

				if(cols.length > col && whatToFind.contains(cols[col])){
					writer.write(line + "\n");
					accepted++;
				}
				if( i % 100000 == 0){
					System.out.println("" + (int)((line.length() * i) / (float)in.length() * 100) + " %,  " + (int)(accepted/(float)(i+1)*100) + " % accepted");
				}												

				i++;
			}

			writer.flush();
			writer.close();
			reader.close();
			System.out.println("DONE in " + (System.currentTimeMillis() - startTime) / 1000 + " seconds");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}
