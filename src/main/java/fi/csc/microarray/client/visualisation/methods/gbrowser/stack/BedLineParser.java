package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

public class BedLineParser extends TsvLineParser {		
	
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
	public BedLineParser(boolean convertCoordinates) {
		this.convertCoordinates = convertCoordinates;
	}

	@Override
	public Region getRegion() {
		
		if (isContentLine()) {

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
	

	
	public String getName() {
		return getString(Column.NAME.ordinal());
	}

	public Float getScore() {
		return getFloat(Column.SCORE.ordinal());
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
	
	public Long getThickStart() {
		return getLong(Column.THICK_START.ordinal());
	}
	
	public Long getThickEnd() {
		return getLong(Column.THICK_END.ordinal());
	}
	
	public String getItemRgb() {
		return getString(Column.ITEM_RGB.ordinal());
	}
	public Long getBlockCount() {
		return getLong(Column.BLOCK_COUNT.ordinal());
	}
	public String getBlockStarts() {
		return getString(Column.BLOCK_STARTS.ordinal());
	}
	
	public int getColumnCount() {
		return values.length;
	}

	@Override
	public String getHeaderStart() {
		return "track";
	}
}
