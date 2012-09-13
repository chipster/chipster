package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Single nucleotide polymorphism (SNP) file parser. File should have a header row.
 * 
 * @author Vilius Zukauskas
 *
 */
public class SNPParser extends TsvParser{

	public SNPParser() {
		super(new FileDefinition(
				Arrays.asList(
						new ColumnDefinition[] {
								new ColumnDefinition(ColumnType.ID, Type.STRING),
								new ColumnDefinition(ColumnType.CHROMOSOME, Type.LONG),
								new ColumnDefinition(ColumnType.POSITION, Type.LONG), //position on chromosome
								new ColumnDefinition(ColumnType.ALLELE, Type.STRING),
								new ColumnDefinition(ColumnType.STRAND, Type.STRING),
								new ColumnDefinition(ColumnType.CONSEQUENCE_TO_TRANSCRIPT, Type.STRING)
						})));
	}
	
	@Override
	public RegionContent[] concise(Chunk chunk) {
		return null;
	}

	@Override
	public String getName() {
		return "SNP Parser";
	}
	
	@Override
	public Region getBpRegion(Chunk chunk) {
		String[] firstRow = getFirstRow(chunk);
	    String[] lastRow = getLastRow(chunk);
		Long start = (Long)get(firstRow, ColumnType.POSITION);
		Long end;
        end = (Long)get(lastRow, ColumnType.POSITION);
        
		Chromosome startChr = (Chromosome)get(firstRow, ColumnType.CHROMOSOME);
		Chromosome endChr = (Chromosome)get(lastRow, ColumnType.CHROMOSOME);
		return new Region(start, startChr, end, endChr);
	}
	
	@Override
	public List<RegionContent> getAll(Chunk chunk,
			Collection<ColumnType> requestedContents) {
		List<RegionContent> rows = new LinkedList<RegionContent>();
		
		for (String row : chunk.getContent().split("\n")) {
			
			LinkedHashMap<ColumnType, Object> values = new LinkedHashMap<ColumnType, Object>();
			
			String[] cols = row.split("\t");
			for (ColumnType requestedContent : requestedContents) {
				
				values.put(requestedContent, this.get(cols, requestedContent));					
			}
			
			Long start = (Long)get(cols, ColumnType.POSITION);
			Long end = (Long)get(cols, ColumnType.POSITION);
			Chromosome chr = (Chromosome)get(cols, ColumnType.CHROMOSOME);
	
			rows.add(new RegionContent(new Region(start, end, chr), values));

		}
		return rows;
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
				return string.equals("-1") || string.equalsIgnoreCase("r") 
				|| string.equals("-") ? Strand.REVERSED	: Strand.FORWARD;

			} else if (col == ColumnType.CHROMOSOME) { 
				return new Chromosome(string);
			} else if (col == ColumnType.CONSEQUENCE_TO_TRANSCRIPT) {
				if (string.startsWith("3PRIME")) {
					return "PRIME3_UTR";
				} else if (string.startsWith("5PRIME")) {
					return "PRIME5_UTR";
				} else if (string.startsWith("STOP_GAINED,FRAMESHIFT")) {
					return "STOP_GAINED_FRAMESHIFT_CODING";
				} else if (string.startsWith("STOP_GAINED,SPLICE_SITE")) {
					return "STOP_GAINED_SPLICED_SITE";
				} else if (string.startsWith("STOP_LOST,SPLICE_SITE")) {
					return "STOP_LOST_SPLICE_SITE";
				} else if (string.startsWith("FRAMESHIFT_CODING,SPLICE_SITE")) {
					return "FRAMESHIFT_CODING_SPLICE_SITE";
				} else if (string.startsWith("STOP_GAINED,FRAMESHIFT_CODING,SPLICE_SITE")) {
					return "STOP_GAINED_FRAMESHIFT_CODING_SPLICE_SITE";
				} else if (string.startsWith("NON_SYNONYMOUS_CODING,SPLICE_SITE")) {
					return "NON_SYNONYMOUS_CODING_SPLICE_SITE";
				} else if (string.startsWith("SPLICE_SITE,SYNONYMOUS_CODING")) {
					return "SPLICE_SITE_SYNONYMOUS_CODING";
				} else if (string.startsWith("SPLICE_SITE,5PRIME_UTR")) {
					return "SPLICE_SITE_5PRIME_UTR";
				} else if (string.startsWith("SPLICE_SITE,3PRIME_UTR")) {
					return "SPLICE_SITE_3PRIME_UTR";
				} else if (string.startsWith("ESSENTIAL_SPLICE_SITE,INTRONIC")) {
					return "ESSENTIAL_SPLICE_SITE_INTRONIC";
				} else if (string.startsWith("SPLICE_SITE,INTRONIC")) {
					return "SPLICE_SITE_INTRONIC";
				} else if (string.startsWith("INTRONIC,NMD_TRANSCRIPT")) {
					return "INTRONIC_NMD_TRANSCRIPT";
				}
				return string;

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
			
		} catch (IndexOutOfBoundsException e) {
			if (col == ColumnType.CONSEQUENCE_TO_TRANSCRIPT) {
				return "NONE";
			} else {
				throw new RuntimeException("error parsing columns: " + Arrays.toString(cols) + " (looking for: " + col + ")", e);
			}
		} catch (Exception e) {
			throw new RuntimeException("error parsing columns: " + Arrays.toString(cols) + " (looking for: " + col + ")", e);
		}
	}

	@Override
	public long getHeaderLength(File file) {
		try {
			FileReader fileReader = new FileReader(file);
			BufferedReader reader = new BufferedReader(fileReader);
			try {
				return reader.readLine().length();
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
		return 0;
	}
}
