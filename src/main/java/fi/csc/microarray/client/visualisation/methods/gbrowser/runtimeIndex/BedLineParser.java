package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.awt.Color;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;

/**
 * This class parses String lines to BedLine objects.
 * 
 * @author klemela
 */
public class BedLineParser extends AbstractTsvLineParser {		
	
	public enum Column {
		
		CHROM ("Chromosome"), 		
		CHROM_START ("Start"), 
		CHROM_END ("End"),
		NAME("Name"), 
		SCORE ("Score"), 
		STRAND ("Strand"), 
		THICK_START ("Thick start"),
		THICK_END ("Thick end"),
		ITEM_RGB ("Item rgb"),
		BLOCK_COUNT ("Block count"),
		BLOCK_SIZES ("Block sizes"),
		BLOCK_STARTS ("Block starts");
		
		private final String name;
		
		Column(String name) {
			this.name = name;
		}
		
		public String toString() {
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

			long start = getLong(Column.CHROM_START.ordinal());
			long end = getLong(Column.CHROM_END.ordinal());
			Chromosome chr = new Chromosome(getString(Column.CHROM.ordinal()));

			if (convertCoordinates) {
				start++;
				end++;
			}

			return new Region(start, end, chr);
			
		} else {
			
			return null;
		}
	}
	
	@Override
	public BedLine getFileLine() {
		BedLine line = new BedLine();
		
		Region region = getRegion();

		line.setChrom(region.start.chr);
		line.setChromStart(region.start.bp);
		line.setChromEnd(region.end.bp);		
		
		int columnCount = values.length;
		int column;
		
		column = Column.NAME.ordinal();
		
		if (columnCount > column) {
			line.setName(getString(column));
		}
		
		column = Column.SCORE.ordinal();
		
		if (columnCount > column) {
			line.setScore(getFloat(column));
		}
		
		column = Column.SCORE.ordinal();
		
		if (columnCount > column) {
			line.setScore(getFloat(column));
		}
		
		column = Column.STRAND.ordinal();
		
		if (columnCount > column) {
			
			String strandString = getString(column);
			
			if ("+".equals(strandString)) {
				line.setStrand(Strand.FORWARD);
			} else if ("-".equals(strandString)) {
				line.setStrand(Strand.REVERSE);
			}
		}
		
		
		column = Column.THICK_START.ordinal();
		
		if (columnCount > column) {
			line.setThickStart(getLong(column));
		}
		
		column = Column.THICK_END.ordinal();
		
		if (columnCount > column) {
			line.setThickEnd(getLong(column));
		}
		
		column = Column.ITEM_RGB.ordinal();
		
		if (columnCount > column) {
			String string = getString(column);
			List<Long> rgb = splitStringToList(string);

			int r = (int)(long)rgb.get(0);
			int g = (int)(long)rgb.get(1);
			int b = (int)(long)rgb.get(2);
			Color c = new Color(r, g, b);
			line.setItemRgb(c);
		}
		
		column = Column.BLOCK_COUNT.ordinal();
		
		if (columnCount > column) {
			line.setBlockCount(getInteger(column));
		}
		
		column = Column.BLOCK_SIZES.ordinal();
		
		if (columnCount > column) {
								
			String string = getString(column);			
			line.setBlockSizes(splitStringToList(string));			
		}
		
		column = Column.BLOCK_STARTS.ordinal();
		
		if (columnCount > column) {
			
			String string = getString(column);					
			line.setBlockStarts(splitStringToList(string));
		}

		return line;
	}
	
	private List<Long> splitStringToList(String string) {
		String[] splitted = string.split(",");
		List<Long> list = new LinkedList<>();
		
		for (String size : splitted) {
			list.add(Long.parseLong(size));
		}
		
		return list;
	}

	public int getColumnCount() {
		return values.length;
	}

	@Override
	public String getHeaderStart() {
		return "track";
	}
}
