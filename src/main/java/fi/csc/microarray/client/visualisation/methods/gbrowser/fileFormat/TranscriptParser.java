package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;

/**
 * Parser for internal trancsript file format.
 * 
 * @author Petri Klemelä
 *
 */
public class TranscriptParser extends AnnotationTsvParser {

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