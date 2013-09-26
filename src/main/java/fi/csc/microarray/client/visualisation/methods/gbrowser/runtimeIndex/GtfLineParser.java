package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.util.HashMap;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;

public class GtfLineParser extends AbstractTsvLineParser {
	
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

	private Map<String, String> attributes;

	@Override
	public Region getRegion() {
		
		if (values == null ){
			return null;
			
		} else {

			long start = getLong(Column.START.ordinal());
			long end = getLong(Column.END.ordinal());
			Chromosome chr = new Chromosome(getString(Column.SEQNAME.ordinal()));
			
			String strandString = getString(Column.STRAND.ordinal());

			Strand strand = null;

			if ("+".equals(strandString)) {
				strand = Strand.FORWARD;
			}

			if ("-".equals(strandString)) {
				strand = Strand.REVERSE;
			}

			return new Region(start, end, chr, strand);
		}
	}

	@Override
	public boolean setLine(String line) {
		
		if (!line.startsWith(getHeaderStart())) {
			this.attributes = null;
		}
		
		return super.setLine(line);
	}

	public String getFeature() {
		return getString(Column.FEATURE.ordinal());
	}
	
	public String getGeneId() {
		return getAttribute("gene_id");
	}
	
	public String getTranscriptId() {
		return getAttribute("transcript_id");
	}
	
	public String getAttribute(String key) {
		if (this.attributes == null) {
			this.attributes = parseAttributes(getString(Column.ATTRIBUTES.ordinal()));			
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

	@Override
	public String getHeaderStart() {
		return "#";
	}

	@Override
	public FileLine getFileLine() {
		// TODO Auto-generated method stub
		return null;
	}
}
