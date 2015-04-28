package fi.csc.microarray.util;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.LinkedHashMap;

import fi.csc.microarray.exception.MicroarrayException;

public class JavaToolUtils {
	
	/**
	 * 
	 * Write additional information about input datasets to a file
	 * 
	 * @param file
	 * @param nameMap
	 * @throws IOException 
	 */
	public static void writeInputDescription(File file, LinkedHashMap<String, String> nameMap) throws IOException {
		
		try (BufferedWriter writer = new BufferedWriter(new FileWriter(file))) {
			writer.write("# Chipster dataset description file\n");
			writer.write("# \n");
			writer.write("# Additional columns may be added later, so don't assume that there will be only two of them.\n");
			writer.write("# Comment lines are allowed only in the beginning of file, but the number of them may vary.\n");
			writer.write("# \n");
			writer.write("# INPUT_NAME	DATASET_NAME\n");
			
			for (String input : nameMap.keySet()) {
				String name = nameMap.get(input);
				writer.write(input + "\t" + name + "\n");
			}
		}
	}
	
	public static LinkedHashMap<String, String> parseOutputDescription(File file) throws IOException, MicroarrayException {
		LinkedHashMap<String, String> nameMap = new LinkedHashMap<>();
		if (file.exists()) {
			try (BufferedReader reader = new BufferedReader(new FileReader(file))) {
				String line;
				while ((line = reader.readLine()) != null) {
					if (line.startsWith("#")) {
						continue;
					}
					String[] splitted = line.split("\t");
					if (splitted.length < 2) {
						throw new MicroarrayException("less than two columns in " + file + " on line '" + line + "'"); 
					}
					String output = splitted[0];
					String name = splitted[1];

					nameMap.put(output, name);
				}
			}
		}
		return nameMap;
	}

	public static void writeOutputDescription(File jobWorkDir, LinkedHashMap<String, String> nameMap) throws IOException {
		/* the format of input and output description files is the same, so 
		 * we can use the same code that writes the input descriptions.
		 */
		writeInputDescription(new File(jobWorkDir, "outputs.tsv"), nameMap);
	}

}
