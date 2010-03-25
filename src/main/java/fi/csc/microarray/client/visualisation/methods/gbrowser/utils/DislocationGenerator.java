package fi.csc.microarray.client.visualisation.methods.gbrowser.utils;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;

public class DislocationGenerator {

	private static final int MAX_BP = 247179968;

	public static void main(String[] args) {

		System.out.println("File generating started");
		Long startTime = System.currentTimeMillis();
		File file = new File("dislocations.fsf");

		BufferedWriter writer;
		try {
			writer = new BufferedWriter(new FileWriter(file));

			Random rand = new Random();

			int i = 0;
			for (int b = 0; b < MAX_BP;) {

				long length = (long) Math.abs((rand.nextGaussian() * 1024));

				String id = "id" + i;
				String fromStart = "" + b;
				String fromEnd = "" + (b + length);
				long toStart = (long) (b + rand.nextGaussian() * MAX_BP / 1028);
				String toEnd = "" + (toStart + length);

				if (b > 0 && toStart > 0 && toStart < MAX_BP) {

					id = fillWithSpaces(id, 16);
					fromStart = fillWithSpaces(fromStart, 16);
					fromEnd = fillWithSpaces(fromEnd, 16);
					String toStartStr = fillWithSpaces("" + toStart, 16);
					String toEndStr = fillWithSpaces("" + toEnd, 50);

					writer.write(id + fromStart + fromEnd + toStartStr + toEndStr + "\n");
				}

				if (i % 100000 == 0) {
					System.out.println("" + (int) ((float) b / MAX_BP * 100) + " % ");
				}

				if (rand.nextFloat() > 0.01) {
					b += Math.pow(rand.nextGaussian(), 2) * 1024;
				} else {
					b += rand.nextGaussian() * 1024 * 1024;
				}

				i++;

			}

			writer.flush();
			writer.close();
			System.out.println("DONE in " + (System.currentTimeMillis() - startTime) / 1000 + " seconds");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static String fillWithSpaces(String orig, int length) {
		String spaces = "                                                                          ";
		return orig + spaces.substring(0, length - orig.length());
	}
}