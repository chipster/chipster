package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Parser for WIG file format.
 * 
 * @see http://genome.ucsc.edu/goldenPath/help/wiggle.html
 */

public class WIGParser extends TsvParser {
	
	private final String VARIABLE_STEP = "variableStep";
	private final String FIXED_STEP = "fixedStep";
	
	
	String type;//variableStep or fixedStep
	String chr;
	Long span;
	Long startPosition;
	Long step;
	
	//fileDefinition is set in the setParser method
	public WIGParser(File file){
		super(null);
		setParser(file);
	}

	@Override
	public RegionContent[] concise(Chunk chunk) {
		return new RegionContent[] {};
	}

	@Override
	public String getName() {
		return "Wiggle file parser";
	}
	
	/**
	 * reading file header info
	 * @param file
	 */
	public void setParser(File file) {
		
		try {
			FileReader fileReader = new FileReader(file);
			BufferedReader reader = new BufferedReader(fileReader);
			String line = reader.readLine();
			
			while (!line.contains("track")){
				line = reader.readLine();
			}
			
			setParser(reader.readLine());
			
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	/**
	 * this chunk is WIG file header info in the middle of data
	 * 
	 * for example
	 *  variableStep chrom=chr1 span=25
	 *   
	 */
	public void setParser(String chunk) {
	
		try {
			String[] cols = chunk.split(" ");
			
			switch (cols.length) {
				case 2:
					type = cols[0];//variable step
					chr = cols[1].replace("chrom=chr", "");
					setFileDefinition(new FileDefinition(Arrays.asList(
							new ColumnDefinition[] { 
									new ColumnDefinition(ColumnType.BP_START, Type.LONG), 
									new ColumnDefinition(ColumnType.VALUE, Type.FLOAT), }
							)));
					break;
				case 3:
					type = cols[0];//variable step
					chr = cols[1].replace("chrom=chr", "");
					span = Long.parseLong(cols[2].replace("span=", ""));
					setFileDefinition(new FileDefinition(Arrays.asList(
							new ColumnDefinition[] { 
									new ColumnDefinition(ColumnType.BP_START, Type.LONG), 
									new ColumnDefinition(ColumnType.VALUE, Type.FLOAT), }
							)));
					break;
				case 4:
					type = cols[0];//fixed step
					chr = cols[1].replace("chrom=chr", "");
					startPosition = Long.parseLong(cols[2].replace("start=", ""));
					step = Long.parseLong(cols[3].replace("step=", ""));
					setFileDefinition(new FileDefinition(Arrays.asList(
							new ColumnDefinition[] { 
									new ColumnDefinition(ColumnType.VALUE, Type.FLOAT), }
							)));
					break;
				case 5:
					type = cols[0];//fixed step
					chr = cols[1].replace("chrom=chr", "");
					startPosition = Long.parseLong(cols[2].replace("start=", ""));
					step = Long.parseLong(cols[3].replace("step=", ""));
					span = Long.parseLong(cols[4].replace("span=", ""));
					setFileDefinition(new FileDefinition(Arrays.asList(
							new ColumnDefinition[] { 
									new ColumnDefinition(ColumnType.VALUE, Type.FLOAT), }
							)));
					break;
				default: 
					break;
			}
		} catch (NumberFormatException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} 
}
	
	@Override
	public List<RegionContent> getAll(Chunk chunk, Collection<ColumnType> requestedContents) {

		List<RegionContent> rows = new LinkedList<RegionContent>();

		if (type.equals(FIXED_STEP)) {
			//fixed step
			
			int nextPosition = 0;
			for (String row : chunk.getContent().split("\n")) {
				Map<ColumnType, Object> values = new HashMap<ColumnType, Object>();
				
				String[] cols = row.split(" ");
				
				if (cols.length > 1) {
					
					setParser(row);
				} else {
					// Calculate start and end positions
					Long start = Long.valueOf(startPosition + nextPosition * step);
					Long end = Long.valueOf(startPosition + nextPosition * step + span - 1);
					
					for (ColumnType requestedContent : requestedContents) {
						values.put(requestedContent, this.get(cols, requestedContent));
					}
					
					rows.add(new RegionContent(new BpCoordRegion(
							start, end, new Chromosome(chr)), values));
				}
				nextPosition++;
			}
			
		} else {
			//variable step
			for (String row : chunk.getContent().split("\n")) {
				
				Map<ColumnType, Object> values = new HashMap<ColumnType, Object>();
				
				String[] cols = row.split("\t");
				
				if (cols.length <2) {
					
					setParser(row);
				} else {

					Long start = Long.parseLong(cols[0]);
					Long end = Long.valueOf(Integer.parseInt(cols[0]) + span-1);
					for (ColumnType requestedContent : requestedContents) {
						values.put(requestedContent, this.get(cols, requestedContent));
					}
			
					rows.add(new RegionContent(new BpCoordRegion(
							start, end, new Chromosome(chr)), values));			
				}
			}	
		}
		
		return rows;
	}
	
	@Override
	public long getHeaderLength(File file) {
		
		FileReader f;
		int lines = 0;
		try {
			
			f = new FileReader(file);
			BufferedReader reader = new BufferedReader(f);
			String str = "";
			String line = "";
			try {
				
				line = reader.readLine();
				while (!isNumber(line.charAt(0))){
					str += line;
					lines++;
					line = reader.readLine();
				}
				
				return str.length() + lines;
				
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		return 0;
	}
	
	private boolean isNumber(char c){
		try {
			Integer.parseInt(String.valueOf(c));
			return true;
			
		} catch (NumberFormatException e) {
			
			return false;
		}
	}
	
	@Override
	public Object get(String[] cols, ColumnType col) {
		
		String string;
		ColumnDefinition fieldDef;
		
		try {

			if (cols.length == 0) {
				return null;
			}
			
			if (col.equals(ColumnType.SEQUENCE)){
				return new String("GGG");
			} else if (col.equals(ColumnType.STRAND)){
				return Strand.BOTH;
			}
			
			string = cols[getFileDefinition().indexOf(col)].trim();
			fieldDef = getFileDefinition().getFieldDef(col);
			
			if (fieldDef.type == Type.FLOAT) {
				return new Float(string);

			} else if (fieldDef.type == Type.LONG) {

				if (string.length() > 0) {
					return new Long(string);
				} else {
					return Long.MIN_VALUE;
				}
			}
			return null;
			
		} catch (Exception e) {
			throw new RuntimeException("error parsing columns: " + Arrays.toString(cols) + " (looking for: " + col + ")", e);
		}
	}
	
	@Override
	public BpCoordRegion getBpRegion(Chunk chunk) {
		
		Long start = 0l;
		Long end = 0l;
		String startChr = chr;
		String endChr;
		String[] header = getLastHeader(chunk);
		
		try {
			endChr = header[0].replace("chrom=chr", "");
		} catch (Exception e) {
			endChr = chr;
		}
						
		if (type.equals(VARIABLE_STEP)) {
			
			try {
				start = Long.valueOf(getFirstRow(chunk)[0]);
				end = Long.valueOf(getLastRow(chunk)[0])+span-1;
			} catch (NumberFormatException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			
		} else {
			
			start = startPosition;
			try {
				end = Long.parseLong(header[1].replace("start=", "")) + 
				getFixedStepEnd(chunk) * step + span - 1;
				
			} catch (Exception e) {
				end = startPosition + getFixedStepEnd(chunk) * step + span - 1;
			}
		}

		return new BpCoordRegion(start, new Chromosome(startChr), end, new Chromosome(endChr));
	}
	
	public String[] getLastHeader(Chunk chunk) {
		
		String content = chunk.getContent();
		
		int lineStartIndex;
		try {
			lineStartIndex = content.lastIndexOf("chrom", content.length() - 2);
		} catch (Exception e) {
			lineStartIndex = -1;
		}
		
		if (lineStartIndex < 0) {
			return null;
		}
		
		return content.substring(lineStartIndex, content.length() - 1).split(" ");
	}
	
	public Long getFixedStepEnd(Chunk chunk) {
		
		int i = 0;
		
		for (String line : chunk.getContent().split("\n")) {
			if (line.contains("chrom")) {
				i=0;
			} else {
				i++;
			}
			
		}
		
		return new Long(i)-1;
	}
	
	@Override
	public String[] getFirstRow(Chunk chunk) {
		
		String content = chunk.getContent();
		
		if (content.startsWith("variableStep")) {
			String row = content.substring(content.indexOf("\n")+1, content.length());
			return content.substring(0, row.indexOf("\n")).split("\t");
			
		} else {
			return content.substring(0, content.indexOf("\n")).split("\t");
			
		}
	}
}
