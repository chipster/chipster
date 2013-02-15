package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.util.HashMap;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

public class StackGtfParser implements Parser {
		
	private static final String GTF_HEADER_START = "#";	
	
	public enum Column {
			
		SEQNAME ("seqname"), 		
		SOURCE ("source"), 
		FEATURE ("feature"),
		START("start"), 
		END ("end"), 
		SCORE ("score"), 
		STRAND ("strand"), 
		FRAME ("frame"), 
		ATTRIBUTES ("attributes");
		
		private final String name;
		
		Column(String name) {
			this.name = name;
		}
		
		String getName() {
			return name;
		}
	}

	private String[] values;
	private Map<String, String> attributes;
	
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
		Chromosome chr = new Chromosome(getString(Column.SEQNAME));
		String strandString = getString(Column.STRAND);
		
		Strand strand = null;
		
		if ("+".equals(strandString)) {
			strand = Strand.FORWARD;
		}
		
		if ("-".equals(strandString)) {
			strand = Strand.REVERSE;
		}
		
		return new Region(start, end, chr, strand);		
	}

	@Override
	public boolean setLine(String line) {
		if (line.startsWith(GTF_HEADER_START)) {
			return false;
		} else {
			this.attributes = null;
			this.values = line.split("\t"); 
			return true; 
		}
	}

	public String getFeature() {
		return getString(Column.FEATURE);
	}
	
	public String getGeneId() {
		return getAttribute("gene_id");
	}
	
	public String getTranscriptId() {
		return getAttribute("transcript_id");
	}
	
	public String getAttribute(String key) {
		if (this.attributes == null) {
			this.attributes = parseAttributes(getString(Column.ATTRIBUTES));			
		}
		return attributes.get(key);
	}

	public static Map<String, String> parseAttributes(String attributeString) {

		String[] split = attributeString.split(";");
		Map<String, String> attributeMap = new HashMap<String, String>(); 

		String key = null;
		String value = null;
		int indexOfQuotationMark = 0;

		for (int i = 0; i < split.length; i++) {

			indexOfQuotationMark = split[i].indexOf("\"");

			key = split[i].substring(0, indexOfQuotationMark - 1).trim();
			value = split[i].substring(indexOfQuotationMark + 1, split[i]
					.lastIndexOf("\""));
			
			attributeMap.put(key, value);
		}

		return attributeMap;
	}
}
