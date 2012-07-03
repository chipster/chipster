package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

/**
 * <p>Parser for BED file format.</p>
 * 
 * <p>BED files describe data lines that are displayed in an annotation track.
 * Information is stored using tab-separated values. Coordinates are converted 
 * from BED's 0-based coordinates to 1-based coordinate system.<br/><br/>
 * 
 * Example:<br/>
 * <pre>
 * track name=pairedReads description="Clone Paired Reads" useScore=1
 * chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
 * chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601
 * </pre></p>
 * 
 * @see http://genome.ucsc.edu/FAQ/FAQformat.html#format1
 * 
 * @author Petri Klemelä, Aleksi Kallio
 *
 */
public class BEDParserWithCoordinateConversion extends BEDParser {

	@Override
	public Object get(String[] cols, ColumnType col) {
		Object obj = super.get(cols, col);
		
		if (col == ColumnType.BP_START || col == ColumnType.BP_END) {
			return (Long)obj + 1;
		}
		return obj;
	}
	
	@Override
	public Region getBpRegion(Chunk chunk) {
		Region reg = super.getBpRegion(chunk);
		return new Region(reg.start.bp + 1, reg.start.chr, reg.end.bp + 1, reg.end.chr);
	}
}