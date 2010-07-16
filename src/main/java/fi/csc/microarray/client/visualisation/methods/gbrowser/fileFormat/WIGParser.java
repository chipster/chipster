package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.LineNumberReader;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Parser for WIG file format.
 * 
 * @see http://genome.ucsc.edu/goldenPath/help/wiggle.html
 */

public class WIGParser extends TsvParser {
	
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
	public RegionContent[] concise(String chunk) {
		return new RegionContent[] {};
	}

	@Override
	public String getName() {
		return "Wiggle file parser";
	}
	
	/**
	 * this chunk is WIG file header info (two lines)
	 * for example 
	 *  browser position chr1:0-1000000
	 *  browser full refGene
	 *  track type=wiggle_0 name="Hs68+FBS PHF8"/ 
	 *  description="H_sapiens-chIPseq-Hs68+FBS PHF8-ALL" /
	 *  visibility=full autoScale=on color=0,200,100 /
	 *  maxHeightPixels=100:50:20 graphType=bar priority=30
	 *  fixedStep chrom=chr1 start=128 step=20 span=20
	 */
	public void setParser(File file) {
		
		try {
			FileReader fileReader = new FileReader(file);
			BufferedReader reader = new BufferedReader(fileReader);
			
			reader.readLine(); //reading browser position chr1:0-1000000 
			reader.readLine(); //reading browser full refGene 
			reader.readLine(); //reading line track type=wiggle_0..
			String[] cols = reader.readLine().split(" ");//splitting the fourth line (variableStep chrom=chr19 ..)
			switch (cols.length) {
				case 2:
					type = cols[0];
					chr = cols[1];
					setFileDefinition(new FileDefinition(Arrays.asList(
							new ColumnDefinition[] { 
									new ColumnDefinition(ColumnType.BP_START, Type.LONG), 
									new ColumnDefinition(ColumnType.VALUE, Type.FLOAT), }
							)));
					break;
				case 3:
					type = cols[0];
					chr = cols[1];
					span = Long.parseLong(cols[2]);
					setFileDefinition(new FileDefinition(Arrays.asList(
							new ColumnDefinition[] { 
									new ColumnDefinition(ColumnType.BP_START, Type.LONG), 
									new ColumnDefinition(ColumnType.VALUE, Type.FLOAT), }
							)));
					break;
				case 4:
					type = cols[0];
					chr = cols[1];
					startPosition = Long.parseLong(cols[2]);
					step = Long.parseLong(cols[3]);
					setFileDefinition(new FileDefinition(Arrays.asList(
							new ColumnDefinition[] { 
									new ColumnDefinition(ColumnType.VALUE, Type.FLOAT), }
							)));
					break;
				case 5:
					type = cols[0];
					chr = cols[1];
					startPosition = Long.parseLong(cols[2]);
					step = Long.parseLong(cols[3]);
					span = Long.parseLong(cols[4]);
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
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	@Override
	public List<RegionContent> getAll(String chunk, Collection<ColumnType> requestedContents) {

		List<RegionContent> rows = new LinkedList<RegionContent>();

		if (type.compareTo("fixedStep") == 0) {
			//fixed step
			for (String row : chunk.split("\n")) {
				Map<ColumnType, Object> values = new HashMap<ColumnType, Object>();
				
				if (requestedContents.contains(ColumnType.CHROMOSOME)) {
					String[] cols = new String [requestedContents.size()];
					
				} else {
					
					
				}
			}
			
		} else {
			//variable step
			for (String row : chunk.split("\n")) {
				
				Map<ColumnType, Object> values = new HashMap<ColumnType, Object>();
				
				String[] cols = row.split("\t");
				
				if (requestedContents.contains(ColumnType.CHROMOSOME)) {
					
					String[] cols2 = new String[3];//two values plus chromosome
					cols2[0] = chr;
					cols2[1] = cols[0];
					cols2[2] = cols[1];
					
					for (ColumnType requestedContent : requestedContents) {
						
						if (span != 1){
							
							Long tmp = Long.parseLong(cols2[0]);
							for (Long i = tmp;i<(tmp+span);i++){
								cols2[1] = i.toString();
								values.put(requestedContent, this.get(cols2, requestedContent));
							}
						} else {
							
							values.put(requestedContent, this.get(cols2, requestedContent));
						}
					}
					
				} else {
					for (ColumnType requestedContent : requestedContents) {
						
						if (span != 1){
							
							Long tmp = Long.parseLong(cols[0]);
							for (Long i = tmp;i<(tmp+span);i++){
								cols[0] = i.toString();
								values.put(requestedContent, this.get(cols, requestedContent));
							}
						} else {
							
							values.put(requestedContent, this.get(cols, requestedContent));
						}
					}
				}
								
				
			}
			
		}
		
		return rows;
	}
	
	@Override
	public long getHeaderLength(File file) {
		
		FileReader f;
		try {
			
			f = new FileReader(file);
			LineNumberReader l = new LineNumberReader(f);
			String str = "";
			try {
				str += l.readLine();
				str += l.readLine();
				
				return str.length() + 2;//two \n
				
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
	
	@Override
	public Object get(String[] cols, ColumnType col) {
		
		try {

			if (cols.length <= 1) {
				return null;
			}

			String string = cols[getFileDefinition().indexOf(col)].trim();

			ColumnDefinition fieldDef = getFileDefinition().getFieldDef(col);

			if (col == ColumnType.STRAND) {
				return string.equalsIgnoreCase("r") || string.equals("-") ? Strand.REVERSED
						: Strand.FORWARD;

			} else if (col == ColumnType.CHROMOSOME) {
				return new Chromosome(string.replace("chr", ""));

			} else if (fieldDef.type == Type.STRING) {
				return string;

			} else if (fieldDef.type == Type.FLOAT) {
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

}
