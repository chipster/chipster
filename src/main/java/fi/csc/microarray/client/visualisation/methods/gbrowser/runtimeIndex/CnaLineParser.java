package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

/**
 * Parser for tsv files produced by CNA tools. 
 * 
 * It is assumed that the header line is given for this parser object before 
 * requesting it to parse any real content rows, because without header line
 * there is no way for parser to know the order and number of columns.  
 * 
 * @author klemela
 */
public class CnaLineParser extends AbstractTsvLineParser {		 
	
	public enum Column {

		ROWNAME ("row name"), 		
		CHROMOSOME ("chromosome"), 
		START ("start"),
		END("end");

		private final String name;

		Column(String name) {
			this.name = name;
		}

		String getName() {
			return name;
		}
	}

	private int gainColumn = -1;
	private int lossColumn = -1;
	private LinkedList<Integer> flagColumns;
	private LinkedList<Integer> logRatioColumns;
	private LinkedList<String> sampleNames;
	
	@Override
	public boolean setLine(String line) {
		if (line.startsWith(getHeaderStart())) {
			parseHeader(line);
			this.values = null;
			return false;
		} else {
			this.values = line.split("\t"); 
			return true; 
		}
	}

	private void parseHeader(String line) {
		
		flagColumns = new LinkedList<Integer>();
		sampleNames = new LinkedList<String>();
		logRatioColumns = new LinkedList<Integer>();
		
		String[] cols = line.split("\t");
		
		for (int i = Column.END.ordinal() + 1; i < cols.length; i++) {
			
			String title = cols[i - 1]; //minus one because header doesn't contain rowname
			
			if (title.startsWith("flag.")) {
				
				this.flagColumns.add(i);				
				
			} else if (title.startsWith("loss.freq")) {
				
				this.lossColumn = i;
				
			} else if (title.startsWith("gain.freq")) {
				
				this.gainColumn = i;
				
			} else if (title.startsWith("segmented.")) {
				
				this.logRatioColumns.add(i);
				this.sampleNames.add(title.replace("segmented.", ""));
			}   			
		}
	}

	@Override
	public Region getRegion() {
		
		if (isContentLine()) {
			
			long start = getLong(Column.START.ordinal());
			long end = getLong(Column.END.ordinal());
			
			Chromosome chr = new Chromosome(getString(Column.CHROMOSOME.ordinal()));
			return new Region(start, end, chr);
			
		} else {
			//This is header line
			return null;
		}
	}
	
	public LinkedList<String> getSampleNames() {
		return sampleNames;
	}
	
	public LinkedList<Float> getFlagValues() {
		
		LinkedList<Float> flagValues = new LinkedList<Float>();
		
		for (int col : flagColumns) {		
			flagValues.add(getFloat(col));
		}
		
		return flagValues;
	}
	
	public LinkedList<Float> getLogRatioValues() {
		
		LinkedList<Float> logRatioValues = new LinkedList<Float>();
		
		for (int col : logRatioColumns) {		
			logRatioValues.add(getFloat(col));
		}
		
		return logRatioValues;
	}
	
	public Float getLossFreq() {
		if (lossColumn != -1) {
			return getFloat(lossColumn);
		}
		return null;
	}
	
	public Float getGainFreq() {
		if (gainColumn != -1) {
			return getFloat(gainColumn);
		}
		return null;
	}
	
	@Override
	public String getHeaderStart() {
		return "chromosome	start	end";
	}

	@Override
	public FileLine getFileLine() {
		// TODO Auto-generated method stub
		return null;
	}
}
