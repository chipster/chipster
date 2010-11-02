package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;

public class ReadPart extends BpCoordRegion {

	
	private String sequencePart;

	public ReadPart() {
		super();
	}

	public ReadPart(BpCoord start, BpCoord end) {
		super(start, end);
	}

	public ReadPart(BpCoordRegion bpRegion) {
		super(bpRegion);
	}

	public ReadPart(Long start, Chromosome chr1, Long end, Chromosome chr2) {
		super(start, chr1, end, chr2);
	}

	public ReadPart(Long start, Long end, Chromosome chr, String sequencePart) {
		super(start, end, chr);
		this.sequencePart = sequencePart;
	}

	public ReadPart(RegionContent read) {
		this(read.region, (String)read.values.get(ColumnType.SEQUENCE));
	}

	public ReadPart(BpCoordRegion region, String sequencePart) {
		this(region);
		this.sequencePart = sequencePart;
	}

	public String getSequencePart() {
		return sequencePart;
	}

	public void setSequencePart(String sequencePart) {
		this.sequencePart = sequencePart;
	}

	
	

}
