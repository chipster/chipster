package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

public class StackBedParser implements Parser {		
		
	private static final String BED_HEADER_START = "track";	
	
	public enum Column {
			
		CHROMOSOME ("chrom"), 		
		START ("chromStart"), 
		END ("chromEnd"),
		NAME("name"), 
		SCORE ("score"), 
		STRAND ("strand"), 
		THICK_START ("thickStart"),
		THICK_END ("thickEnd"),
		ITEM_RGB ("itemRgb"),
		BLOCK_COUNT ("blockCount"),
		BLOCK_SIZES ("blockSizes"),
		BLOCK_STARTS ("blockStarts");
		
		private final String name;
		
		Column(String name) {
			this.name = name;
		}
		
		String getName() {
			return name;
		}
	}

	private String[] values;
	private boolean convertCoordinates;
	
	/** 
	 * @param convertCoordinates Convert coordinates from bed file 0-based coordinate system to genome browser 1-based coordinates.
	 */
	public StackBedParser(boolean convertCoordinates) {
		this.convertCoordinates = convertCoordinates;
	}
	
	private long getLong(Column column) {
		String string = values[column.ordinal()];
		return new Long(string);
	}
	
	private String getString(Column column) {
		return values[column.ordinal()];
	}

	@Override
	public Region getRegion() {
			
		long start = getLong(Column.START);
		long end = getLong(Column.END);
		Chromosome chr = new Chromosome(getString(Column.CHROMOSOME));
		
		if (convertCoordinates) {
			start++;
			end++;
		}
		
		return new Region(start, end, chr);		
	}

	@Override
	public boolean setLine(String line) {
		if (line.startsWith(BED_HEADER_START)) {
			return false;
		} else {
			this.values = line.split("\t"); 
			return true; 
		}
	}
	
	public Strand getStrand() {
		String strandString = getString(Column.STRAND);
		
		Strand strand = null;
		
		if ("+".equals(strandString)) {
			strand = Strand.FORWARD;
		}
		
		if ("-".equals(strandString)) {
			strand = Strand.REVERSE;
		}
		
		return strand;
	}

	public String getScore() {
		return getString(Column.SCORE);
	}
}
