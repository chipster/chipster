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

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
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
	 *  track type=wiggle_0 name="H9 Oct4 replicate 2" description="H9 Oct4/
	 *   replicate 2 in reads per million" visibility=full autoScale=on color=0,0,0
	 *   variableStep chrom=chr1 span=25
	 */
	public void setParser(File file) {
		
		try {
			FileReader fileReader = new FileReader(file);
			BufferedReader reader = new BufferedReader(fileReader);
			String[] cols = null;
			//TODO change
//			reader.readLine(); //reading browser position chr1:0-1000000 
//			reader.readLine(); //reading browser full refGene 
//			reader.readLine(); //reading line track type=wiggle_0..
			String line = reader.readLine();
			
			if (line.contains("track")){
				cols = reader.readLine().split(" ");//splitting line (variableStep chrom=chr1 ..)
				
			} else if (line.contains("chrom=")){
				//splitting line 
				cols = line.split(" ");
			}
			
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
		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

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
	public List<RegionContent> getAll(String chunk, Collection<ColumnType> requestedContents) {

		List<RegionContent> rows = new LinkedList<RegionContent>();

		if (type.compareTo("fixedStep") == 0) {
			//fixed step
			
			int nextPosition = 0;
			for (String row : chunk.split("\n")) {
				Map<ColumnType, Object> values = new HashMap<ColumnType, Object>();
				
				String[] cols = row.split(" ");
				
				if (cols.length > 1) {
					
					setParser(row);
				} else {
					
					//bp_start,bp_end,chromosome,value
					String[] cols2 = new String[4];
					cols2[0] = chr;
					cols2[1] = String.valueOf(startPosition + nextPosition*step); //bp_start
					cols2[2] = row;//value
					
					//bp_end
					if (span != 1) {
						cols2[3] = String.valueOf(startPosition + nextPosition*step
								+ span-1);
					} else {
						cols2[3] = String.valueOf(startPosition + nextPosition*step);
					}
					
					for (ColumnType requestedContent : requestedContents) {
						
						values.put(requestedContent, this.get(cols2, requestedContent));
					}
					
					Long start = (Long)get(cols2, ColumnType.BP_START);
					Long end = (Long)get(cols2, ColumnType.BP_END);
					Chromosome chr = (Chromosome)get(cols2, ColumnType.CHROMOSOME);
			
					rows.add(new RegionContent(new BpCoordRegion(start, end, chr), values));
				}
				nextPosition++;
			}
			
		} else {
			//variable step
			for (String row : chunk.split("\n")) {
				
				Map<ColumnType, Object> values = new HashMap<ColumnType, Object>();
				
				String[] cols = row.split("\t");
				if (cols.length <2) {
					
					setParser(row);
				} else {
					
					//two values(bp_start, value) plus chromosome, and bp_end
					String[] cols2 = new String[4];
					cols2[0] = chr;
					cols2[1] = cols[0];
					cols2[2] = cols[1];
					
					if (span != 1) {
						cols2[3] = String.valueOf(Integer.parseInt(cols[0]) + span-1);
					} else {
						cols2[3] = cols[0];
					}
					
					for (ColumnType requestedContent : requestedContents) {
						
						values.put(requestedContent, this.get(cols2, requestedContent));
					}
					
					Long start = (Long)get(cols2, ColumnType.BP_START);
					Long end = (Long)get(cols2, ColumnType.BP_END);
					Chromosome chr = (Chromosome)get(cols2, ColumnType.CHROMOSOME);
			
					rows.add(new RegionContent(new BpCoordRegion(start, end, chr), values));			
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
				//TODO change
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
		
		String string;
		ColumnDefinition fieldDef;
		
		try {

			if (cols.length <= 1) {
				return null;
			}

			if (col == ColumnType.CHROMOSOME) {
				return new Chromosome(cols[0].replace("chr", ""));

			} else if (col == ColumnType.BP_END){
				//System.out.println(cols[3]);
				return new Long(cols[3]);
				
			} else if (col == ColumnType.BP_START){
				return new Long(cols[1]);
				
			} else {
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
			}
			return null;
			
		} catch (Exception e) {
			throw new RuntimeException("error parsing columns: " + Arrays.toString(cols) + " (looking for: " + col + ")", e);
		}
		
	}

}
