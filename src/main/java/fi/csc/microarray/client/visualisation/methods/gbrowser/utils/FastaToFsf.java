package fi.csc.microarray.client.visualisation.methods.gbrowser.utils;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

public class FastaToFsf {

	public static void main(String[] args) {

		System.out.println("File converting started");
		Long startTime = System.currentTimeMillis();

		File in = new File("chr1.fa");
		File out = new File("chr1.fsf");

		BufferedWriter writer;
		BufferedReader reader;
		try {
			writer = new BufferedWriter(new FileWriter(out));
			reader = new BufferedReader(new FileReader(in));
			// Header
			reader.readLine();

			int i = 0;

			long bp = 1;
			String line;

			while ((line = reader.readLine()) != null) {

				StringBuilder constLine = new StringBuilder();

				constLine.append(fillWithSpaces(line, 64));

				constLine.append(fillWithSpaces("" + bp, 16));

				bp += line.length();

				constLine.append("\n");

				writer.write(constLine.toString());

				if (i % 10000 == 0) {
					System.out.println("" + (int) ((line.length() * i) / (float) in.length() * 100) + " % ");
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

	private static String fillWithSpaces(String orig, int length) {
		String spaces = "                                                                                 ";
		return orig + spaces.substring(0, length - orig.length());
	}
}
