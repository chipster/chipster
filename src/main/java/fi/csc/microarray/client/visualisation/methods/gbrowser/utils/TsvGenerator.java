package fi.csc.microarray.client.visualisation.methods.gbrowser.utils;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Random;


public class TsvGenerator {	
	

	private static long rows = (long)Math.pow(10, 7);
	
	public static void main(String[] args) {
		
		System.out.println("File generating started");
		Long start = System.currentTimeMillis();
		File file = new File("tsv.tsv");		

		BufferedWriter writer;
		try {
			writer = new BufferedWriter(new FileWriter(file));

			Random rand =  new Random();
			
			writer.write(" " + "\t" + "chip.microarray001" + "\t" + "annotation" + "\t" + "description" + "\n");

			
			for (int i = 0; i < rows; i++){

				String id = "" + i;
				String value = "" + calcY(i);
				String strand = rand.nextBoolean()? "+" : "-";
				String seq = "";

				for (int j = 0; j < 16; j++){
					int r = rand.nextInt(4);
					switch(r){
					case 0: seq += "A"; break;
					case 1: seq += "T"; break;
					case 2: seq += "G"; break;
					case 3: seq += "C"; break;
					}
				}

//				id = fillWithSpaces(id, 16);
//				value = fillWithSpaces(value, 16);
//				strand = fillWithSpaces(strand, 2);

				writer.write(id + "\t" + value + "\t" + strand + "\t" + seq + "\n");
				if( i % 10000 == 0){
					System.out.println(i);
				}
			}
			
			writer.flush();
			writer.close();
			System.out.println("DONE in " + (System.currentTimeMillis() - start) / 1000 + " seconds");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private static String fillWithSpaces(String orig, int length){
		String spaces = "                            ";
		return orig + spaces.substring(0, length - orig.length());
	}

	private static float calcY(int x){
		float value = 0;
		
		for(float factor = rows / 3; factor > 3; factor /= 2){
			value += Math.sin(x / (factor)) * Math.log(factor);			
		}		
		
		return value;
	}
}