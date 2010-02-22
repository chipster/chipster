package fi.csc.microarray.client.visualisation.methods.gbrowser.utils;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;


public class IntervalFileGenerator {
	
	private static final int MAX_BP = 247179968;

	public static void main(String[] args) {
		
		System.out.println("File generating started");
		Long startTime = System.currentTimeMillis();
		File file = new File("reads.fsf");		

		BufferedWriter writer;
		try {
			writer = new BufferedWriter(new FileWriter(file));

			Random rand =  new Random();

			int i = 0;
			for (int b = 0; b < MAX_BP; ){

				String id = "id" + i;
				String strand = rand.nextBoolean()? "+" : "-";
				String ref = "" + rand.nextInt(100);
				String start = "" + b;
				
				int length = rand.nextInt(25) + 25;
				StringBuilder seq = new StringBuilder();

				for (int j = 0; j < length; j++){
					seq.append(getRandomChar(rand));
				}
				
				StringBuilder quality = new StringBuilder();
				for(int j = 0; j < length; j++){
					quality.append("I");
				}
				
				String other = "" + rand.nextInt(1000);
				String mismatch = "" + rand.nextInt(50) + getRandomChar(rand) + ">" + getRandomChar(rand);
				

				id = fillWithSpaces(id, 16);
				strand = fillWithSpaces(strand, 2);
				ref = fillWithSpaces(ref, 16);
				start = fillWithSpaces(start, 16);			
				String seqStr = fillWithSpaces(seq.toString(), 50);
				String qualStr = fillWithSpaces(quality.toString(), 50);
				other = fillWithSpaces(other, 16);
				mismatch = fillWithSpaces(mismatch, 16);				

				writer.write(id + strand + ref + start + seqStr + qualStr + other + mismatch + "\n");
				if( i % 100000 == 0){
					System.out.println("" + (int)((float)b / MAX_BP * 100) + " % ");
				}
				
				if(rand.nextFloat() > 0.0001){
					b += Math.pow(rand.nextGaussian(), 2) * 50;
				} else {
					b += rand.nextInt(1024*1024);
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

	private static String getRandomChar(Random rand) {
		int r = rand.nextInt(4);
		switch(r){
		case 0: return "A"; 
		case 1: return "T"; 
		case 2: return "G"; 
		case 3: return "C"; 
		}
		return null;
	}

	private static String fillWithSpaces(String orig, int length){
		String spaces = "                                                                          ";
		return orig + spaces.substring(0, length - orig.length());
	}
}