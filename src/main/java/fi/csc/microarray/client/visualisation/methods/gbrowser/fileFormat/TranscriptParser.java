package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class TranscriptParser extends ConcisedTsvParser {

	public TranscriptParser() {
		super(new FileDefinition(
				Arrays.asList(
						new ColumnDefinition[] {
								new ColumnDefinition(ColumnType.CHROMOSOME, Type.LONG),
								new ColumnDefinition(ColumnType.PARENT_BP_START, Type.LONG),
								new ColumnDefinition(ColumnType.PARENT_BP_END, Type.LONG),
								new ColumnDefinition(ColumnType.STRAND, Type.STRING),
								new ColumnDefinition(ColumnType.BP_START, Type.LONG),
								new ColumnDefinition(ColumnType.BP_END, Type.LONG),
								new ColumnDefinition(ColumnType.DESCRIPTION, Type.STRING),
								new ColumnDefinition(ColumnType.VALUE, Type.STRING),							
								new ColumnDefinition(ColumnType.PARENT_ID, Type.STRING),
								new ColumnDefinition(ColumnType.PARENT_PART, Type.STRING)
		})));
	}

	@Override
	public String getName() {
		return "Chipster transcript file";
	}
	
	@Override
	public long getDefaulChunkLength() {
		return 512;
	}
}