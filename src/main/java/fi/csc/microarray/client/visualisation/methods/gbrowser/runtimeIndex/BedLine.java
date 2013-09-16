package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.awt.Color;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.BedLineParser.Column;

/**
 * A class which represents single line of bed file.
 * 
 * @author klemela
 */
public class BedLine extends FileLine implements ScatterplotValue, SelectionText  {
	
	private Chromosome chrom;
	private Long chromStart;
	private Long chromEnd;
	private String name;
	private Float score;
	private Strand strand;
	private Long thickStart;
	private Long thickEnd;
	private Color itemRgb;
	private Integer blockCount;
	private List<Long> blockSizes;
	private List<Long> blockStarts;	
	
	public Chromosome getChrom() {
		return chrom;
	}
	public void setChrom(Chromosome chromosome) {
		this.chrom = chromosome;
	}
	public Long getChromStart() {
		return chromStart;
	}
	public void setChromStart(Long chromStart) {
		this.chromStart = chromStart;
	}
	public Long getChromEnd() {
		return chromEnd;
	}
	public void setChromEnd(Long chromEnd) {
		this.chromEnd = chromEnd;
	}
	public String getName() {
		return name;
	}
	public void setName(String name) {
		this.name = name;
	}
	public Float getScore() {
		return score;
	}
	public void setScore(Float score) {
		this.score = score;
	}
	public Strand getStrand() {
		return strand;
	}
	public void setStrand(Strand strand) {
		this.strand = strand;
	}
	public Long getThickStart() {
		return thickStart;
	}
	public void setThickStart(Long thickStart) {
		this.thickStart = thickStart;
	}
	public Long getThickEnd() {
		return thickEnd;
	}
	public void setThickEnd(Long thickEnd) {
		this.thickEnd = thickEnd;
	}
	public Color getItemRgb() {
		return itemRgb;
	}
	public void setItemRgb(Color itemRgb) {
		this.itemRgb = itemRgb;
	}
	public Integer getBlockCount() {
		return blockCount;
	}
	public void setBlockCount(Integer blockCount) {
		this.blockCount = blockCount;
	}
	public List<Long> getBlockSizes() {
		return blockSizes;
	}
	public void setBlockSizes(List<Long> blockSizes) {
		this.blockSizes = blockSizes;
	}
	public List<Long> getBlockStarts() {
		return blockStarts;
	}
	public void setBlockStarts(List<Long> blockStarts) {
		this.blockStarts = blockStarts;
	}
	
	@Override
	public Float getScatterplotValue() {
		return getScore();
	}
	@Override
	public Color getScatterplotColor() {
		return getItemRgb();
	}
	@Override
	public String getText() {			
		
		String[] columns = new String[] {
				"" + Column.CHROM, "" + getChrom(),
				"" + Column.CHROM_START, "" + getChromStart(), 
				"" + Column.CHROM_END, "" + getChromEnd(),
				"" + Column.NAME, "" + getName(),
				"" + Column.SCORE, ""  + getScore(),
				"" + Column.STRAND, "" + getStrand(),
				"" + Column.THICK_START, "" + getThickStart(),
				"" + Column.THICK_END, "" + getThickEnd(),
				"" + Column.ITEM_RGB, "" + getItemRgbAsString(),
				"" + Column.BLOCK_COUNT, "" + getBlockCount(),
				"" + Column.BLOCK_SIZES, "" + getBlockSizes(),
				"" + Column.BLOCK_STARTS, "" + getBlockStarts()				
		};
		
		return super.format("BED feature info", columns);
	}
	
	private String getItemRgbAsString() {
		Color color = getItemRgb();
		String colorText = null;
		if (color != null) {
			colorText = "r="+color.getRed() + ", g=" + color.getGreen() + ", b=" + color.getBlue();
		}
		return colorText;
	}	
}
