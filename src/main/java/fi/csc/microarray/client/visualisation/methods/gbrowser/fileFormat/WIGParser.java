package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Parser for the WIG file format.
 * 
 * @author Vilius Zukauskas
 * @see http://genome.ucsc.edu/goldenPath/help/wiggle.html
 */
public class WIGParser extends TsvParser {

	public WIGParser() {
		super(new FileDefinition(
				Arrays.asList(
						new ColumnDefinition[] {
								new ColumnDefinition(ColumnType.CHROMOSOME, Type.LONG),
								new ColumnDefinition(ColumnType.BP_START, Type.LONG),
								new ColumnDefinition(ColumnType.BP_END, Type.LONG),
								new ColumnDefinition(ColumnType.VALUE, Type.STRING)
						})));
	}

	public WIGParser(FileDefinition fileDefinition) {
		super(fileDefinition);
	}   

	@Override
	public RegionContent[] concise(Chunk chunk) {
		return null;
	}

	@Override
	public String getName() {
		return "Wiggle parser";
	}
	
	
}
