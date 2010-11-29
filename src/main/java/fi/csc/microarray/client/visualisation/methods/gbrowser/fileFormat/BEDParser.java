package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Parser for BED file format.
 * 
 * BED files describe data lines that are displayed in an annotation track.
 * Information is stored using tab-separated values.<br/><br/>
 * 
 * Example:<br/>
 * <pre>
 * track name=pairedReads description="Clone Paired Reads" useScore=1
 * chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512
 * chr22 2000 6000 cloneB 900 - 2000 6000 0 2 433,399, 0,3601
 * </pre>
 * 
 * @see http://genome.ucsc.edu/FAQ/FAQformat.html#format1
 *
 */
public class BEDParser extends TsvParser {

	public BEDParser() {
		super(new FileDefinition(Arrays.asList(
				new ColumnDefinition[] { 
						new ColumnDefinition(ColumnType.CHROMOSOME, Type.STRING), 
						new ColumnDefinition(ColumnType.BP_START, Type.LONG), 
						new ColumnDefinition(ColumnType.BP_END, Type.LONG), }
				)));
	}

	public BEDParser(FileDefinition fileDefinition) {
		super(fileDefinition);
	}

	@Override
	public String getName() {
		return "Chipster peaks";
	}

	@Override
	public long getDefaulChunkLength() {
		return 128;
	}

	@Override
	public RegionContent[] concise(Chunk chunk) {
		return new RegionContent[] {};
	}

	@Override
	public Object get(String[] cols, ColumnType col) {

		Object obj = super.get(cols, col);
		
		return obj;
	}
}