package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

public class CytobandLineParser extends AbstractTsvLineParser {		 
	
	public enum Column {

		CHR ("chromosome"), 		
		KARYOTYPE_ID ("karyotype id"), 
		START ("start"),
		END("end"), 
		BAND ("band"), 
		STAIN ("stain");

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
				
		long start = getLong(Column.START.ordinal());
		long end = getLong(Column.END.ordinal());
			
		Chromosome chr = new Chromosome(getString(Column.CHR.ordinal()));
		
		return new Region(start, end, chr);		
	}

	public String getBand() {
		return getString(Column.BAND.ordinal());
	}
	
	public String getStain() {
		return getString(Column.STAIN.ordinal());
	}

	@Override
	public FileLine getFileLine() {
		// TODO Auto-generated method stub
		return null;
	}
}
