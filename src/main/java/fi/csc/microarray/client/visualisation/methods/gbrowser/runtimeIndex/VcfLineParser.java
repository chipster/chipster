package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

public class VcfLineParser extends AbstractTsvLineParser {		 
	
	public enum Column {

		CHROM ("CHROM"), 		
		POS ("POS"), 
		ID ("ID"),
		REF("REF"), 
		ALT ("ALT"), 
		QUAL ("QUAL"), 
		FILTER ("FILTER"),
		INFO ("INFO");

		private final String name;

		Column(String name) {
			this.name = name;
		}

		public String toString() {
			return name;
		}
	}

	@Override
	public Region getRegion() {
		
		if (isContentLine()) {
			
			long start = getLong(Column.POS.ordinal());
			
			Chromosome chr = new Chromosome(getString(Column.CHROM.ordinal()));
			return new Region(start, start, chr);
			
		} else {
			//This is header line
			return null;
		}
	}

	public VcfLine getFileLine() {
		VcfLine line = new VcfLine();
		
		Region region = getRegion();
		
		line.setChrom(region.start.chr);
		line.setPos(region.start.bp);
		line.setId(getString(Column.ID.ordinal()));
		line.setRef(getString(Column.REF.ordinal()));
		line.setAlt(getString(Column.ALT.ordinal()));
		line.setQual(getFloat(Column.QUAL.ordinal()));
		line.setFilter(getString(Column.FILTER.ordinal()));
		line.setInfo(getString(Column.INFO.ordinal()));
		
		return line;
	}

	@Override
	public String getHeaderStart() {
		return "#";
	}
}
