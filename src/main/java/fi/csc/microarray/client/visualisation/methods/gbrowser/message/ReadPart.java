package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;

/**
 * Single spliced part of a spliced read. 

 * @author Aleksi Kallio
 *
 */
public class ReadPart extends Region {
	
	private String sequencePart;
	private RegionContent read;

	public ReadPart(Long start, Long end, Chromosome chr, RegionContent read, String sequencePart) {
		super(start, end, chr);
		this.read = read;
		this.sequencePart = sequencePart;
	}

	public ReadPart(RegionContent read) {
		super(read.region);
		this.read = read;
		this.sequencePart = (String)read.values.get(ColumnType.SEQUENCE);
	}

	public RegionContent getRead() {
		return read;
	}

	public String getSequencePart() {
		return sequencePart;
	}

	public void setSequencePart(String sequencePart) {
		this.sequencePart = sequencePart;
	}

	
	

}
