package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

public class StackBedParser extends StackTsvParser {		
	
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
	
	private boolean convertCoordinates;
	
	/** 
	 * @param convertCoordinates Convert coordinates from bed file 0-based coordinate system to genome browser 1-based coordinates.
	 */
	public StackBedParser(boolean convertCoordinates) {
		this.convertCoordinates = convertCoordinates;
	}

	@Override
	public Region getRegion() {
		
		if (isContentLinle()) {

			long start = getLong(Column.START.ordinal());
			long end = getLong(Column.END.ordinal());
			Chromosome chr = new Chromosome(getString(Column.CHROMOSOME.ordinal()));

			if (convertCoordinates) {
				start++;
				end++;
			}

			return new Region(start, end, chr);
			
		} else {
			
			return null;
		}
	}
	
	public Strand getStrand() {
		String strandString = getString(Column.STRAND.ordinal());
		
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
		return getString(Column.SCORE.ordinal());
	}

	@Override
	public String getHeaderStart() {
		return "track";
	}
}
