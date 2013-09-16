package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.VcfLineParser.Column;


/**
 * This class represents a single line of Vcf file.
 * 
 * @author klemela
 */
public class VcfLine extends FileLine implements SelectionText  {

	private Chromosome chrom;
	private Long pos;
	private String id;
	private String ref;
	private String alt;
	private Float qual;
	private String filter;
	private String info;	

	public Chromosome getChrom() {
		return chrom;
	}
	public void setChrom(Chromosome chrom) {
		this.chrom = chrom;
	}
	public Long getPos() {
		return pos;
	}
	public void setPos(Long pos) {
		this.pos = pos;
	}
	public String getId() {
		return id;
	}
	public void setId(String id) {
		this.id = id;
	}
	public String getRef() {
		return ref;
	}
	public void setRef(String ref) {
		this.ref = ref;
	}
	public String getAlt() {
		return alt;
	}
	public void setAlt(String alt) {
		this.alt = alt;
	}
	public Float getQual() {
		return qual;
	}
	public void setQual(Float qual) {
		this.qual = qual;
	}
	public String getFilter() {
		return filter;
	}
	public void setFilter(String filter) {
		this.filter = filter;
	}
	public String getInfo() {
		return info;
	}
	public void setInfo(String info) {
		this.info = info;
	}	

	@Override
	public String getText() {				

		String[] columns = new String[] {
				"" + Column.CHROM, "" + getChrom(),
				"" + Column.POS, "" + getPos(), 
				"" + Column.ID, "" + getId(),
				"" + Column.REF, "" + getRef(),
				"" + Column.ALT, ""  + getAlt(),
				"" + Column.QUAL, "" + getQual(),
				"" + Column.FILTER, "" + getFilter(),
				"" + Column.INFO, "" + formatList(getInfo()),				
		};
		
		return super.format("VCF variant info", columns);
	}
	
	private String formatList(String list) {
		StringBuilder builder = new StringBuilder();
		builder.append("\n");
		for (String item : list.split(";")) {
			builder.append("\t");
			builder.append(item);
			builder.append("\n");
		}
		return builder.toString();
	}
}
