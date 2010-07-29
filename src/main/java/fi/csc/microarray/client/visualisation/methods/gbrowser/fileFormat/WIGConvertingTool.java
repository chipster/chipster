package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;

public class WIGConvertingTool {
	
	private class HeaderDefinition {
		
		String type;//variableStep or fixedStep
		String chr;
		Long span;
		Long startPosition;
		Long step;
		Long headerPosition;
	}
	
	HeaderDefinition header;
	
	public void convert(File file){
		
		String line;
		int linesCount = 0;
		
		FileReader fReader;
		try {
			fReader = new FileReader(file);
			BufferedReader reader = new BufferedReader(fReader);
			
			while (reader.ready()) {
				line = reader.readLine();
				if (line.contains("Step")){
					header = setHeader(line);
				} else {
					line = formLine(line, linesCount);
					//write line
				}
			}
			
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
					header.type = cols[0];//variable step
					header.chr = cols[1].replace("chrom=chr", "");
					
					break;
				case 3:
					header.type = cols[0];//variable step
					header.chr = cols[1].replace("chrom=chr", "");
					header.span = Long.parseLong(cols[2].replace("span=", ""));
					break;
				case 4:
					header.type = cols[0];//fixed step
					header.chr = cols[1];//.replace("chrom=chr", "");
					header.startPosition = Long.parseLong(cols[2].replace("start=", ""));
					header.step = Long.parseLong(cols[3].replace("step=", ""));
					break;
				case 5:
					header.type = cols[0];//fixed step
					header.chr = cols[1];//.replace("chrom=chr", "");
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
		//TODO formLine
	}

}
