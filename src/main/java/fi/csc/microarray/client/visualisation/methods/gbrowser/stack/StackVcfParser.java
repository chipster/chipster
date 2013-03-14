package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

public class StackVcfParser extends StackTsvParser {		 
	
	public enum Column {

		CHROM ("chromosome"), 		
		POS ("position"), 
		ID ("identifiers"),
		REF("references base(s)"), 
		ALT ("non-reference alleles"), 
		QUAL ("quality score"), 
		FILTER ("filter pass"),
		INFO ("additional information");

		private final String name;

		Column(String name) {
			this.name = name;
		}

		String getName() {
			return name;
		}
	}

	@Override
	public Region getRegion() {
		
		if (isContentLinle()) {
			
			long start = getLong(Column.POS.ordinal());
			
			Chromosome chr = new Chromosome(getString(Column.CHROM.ordinal()));
			return new Region(start, start, chr);
			
		} else {
			//This is header line
			return null;
		}
	}

	public String getQuality() {
		return getString(Column.QUAL.ordinal());
	}

	@Override
	public String getHeaderStart() {
		return "#";
	}
}
