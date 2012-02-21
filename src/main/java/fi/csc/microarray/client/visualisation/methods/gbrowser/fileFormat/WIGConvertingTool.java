package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 * 
 * Simple wig format conversion tool.
 * 
 * @author Vilius Zukauskas
 *  
 */
public class WIGConvertingTool {
	
	private static final String VARIABLE_STEP = "variableStep";
	private static final String FIXED_STEP = "fixedStep";
	
	private class HeaderDefinition {
		
		String type;//variableStep or fixedStep
		String chr;
		Long span;
		Long startPosition;
		Long step;
	}
	
	HeaderDefinition header;
	float value = -1;
	String line = "";
	String previousLine = "";
	String output;
	int linesCount = 1;
	
	public void convert(File file) {
				
		FileReader fReader;
		BufferedWriter writer;
		try {
			fReader = new FileReader(file);
			BufferedReader reader = new BufferedReader(fReader);
			
			writer = new BufferedWriter(new FileWriter(file.getAbsolutePath()+".out"));
			
			while (!line.contains("track")) {
				line = reader.readLine();
			}
			
			//read the first header, and the first value line
			
			line = reader.readLine();
			header = setHeader(line);
			line = reader.readLine();
			
			if (header.type.equals(VARIABLE_STEP)) {
                header.startPosition = Long.parseLong(line.split("\t")[0]);
			}

			
			value = getValue(line);
			
			//TODO optimize
			while (reader.ready()) {
				previousLine = line;
				line = reader.readLine();
				if (line.contains("Step")) {
					
					if (linesCount != 0) {
						output = formLine(previousLine, linesCount-1);
						writer.write(output);
					}
					header = setHeader(line);
					linesCount = 0;
					
					previousLine = line;
					line = reader.readLine();
					if (header.type.equals(VARIABLE_STEP)) {
						header.startPosition = Long.parseLong(line.split("\t")[0]);
						value = getValue(line);
						
					} else {
						value = Long.parseLong(line);
						
					}
					
				} else {
					if (header.type.equals(FIXED_STEP)) {
						if (value != getValue(line)) {
							output = formLine(previousLine, linesCount);
							writer.write(output);
							header.startPosition = header.startPosition + linesCount * header.step;
							linesCount = 0;
							value = getValue(line);
						}
					} else {
						
						if (getVarStepPos(line) != (getVarStepPos(previousLine) + header.span)) { 
							output = formLine(previousLine, linesCount);
							writer.write(output);
							linesCount = 0;
							value = getValue(line);
							header.startPosition = getVarStepPos(line);
						} else {
							if (value != getValue(line)) {
								output = formLine(previousLine, linesCount);
								writer.write(output);
								linesCount = 0;
								value = getValue(line);
								header.startPosition = getVarStepPos(line);
							}
						}
					}
					linesCount++;
				}
			}

			output = formLine(line, linesCount);
			writer.write(output);
			
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	public HeaderDefinition setHeader(String chunk) {

		HeaderDefinition header = new HeaderDefinition();
		
		try {
			
			if (chunk.indexOf("\n") != -1){
				chunk = chunk.substring(0, chunk.indexOf("\n"));
			}
			String[] cols = chunk.split(" ");
			
			//TODO configure depending on converting
			switch (cols.length) {
				case 2:
					header.type = cols[0]; //variable step
					header.chr = cols[1].replace("chrom=chr", "");
					
					break;
				case 3:
					header.type = cols[0]; //variable step
					header.chr = cols[1].replace("chrom=chr", "");
					header.span = Long.parseLong(cols[2].replace("span=", ""));
					break;
				case 4:
					header.type = cols[0]; //fixed step
					header.chr = cols[1].replace("chrom=", "");
					header.startPosition = Long.parseLong(cols[2].replace("start=", ""));
					header.step = Long.parseLong(cols[3].replace("step=", ""));
					break;
				case 5:
					header.type = cols[0]; //fixed step
					header.chr = cols[1].replace("chrom=", "");
					header.startPosition = Long.parseLong(cols[2].replace("start=", ""));
					header.step = Long.parseLong(cols[3].replace("step=", ""));
					header.span = Long.parseLong(cols[4].replace("span=", ""));
					break;
				default: 
					break;
			}
			
			return header;
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			return null;
		}
	}
	
	public String formLine(String line, int linesCount) {
		
		if (header.type.startsWith("variable")) {
			String[] splitted = line.split("\t");
			
			return header.chr + "\t" + header.startPosition + "\t" + 
			(Long.parseLong(splitted[0]) + header.span - 1) + "\t" + splitted[1] + "\n";  
		} else {
			return header.chr + "\t" + header.startPosition + "\t" + 
			(header.startPosition + linesCount * header.span - 1) + "\t" + line + "\n";
		}
	}
	
	public float getValue(String value){
		String[] splitted = value.split("\t"); 
		switch (splitted.length) {
		case 1:
			return Float.parseFloat(splitted[0]);
		case 2:
			return Float.parseFloat(splitted[1]);
		default:
			return (Float)null;
		}
	}
	
	public long getVarStepPos(String value) {
		return Long.parseLong(value.split("\t")[0]);
	}

}
