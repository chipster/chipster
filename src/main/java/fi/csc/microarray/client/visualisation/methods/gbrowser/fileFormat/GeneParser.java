package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;

public class GeneParser extends ConcisedTsvParser{

	public GeneParser() {
		super(new FileDefinition(
				Arrays.asList(
						new ColumnDefinition[] {
								new ColumnDefinition(ColumnType.CHROMOSOME, Type.LONG),
								new ColumnDefinition(ColumnType.BP_START, Type.LONG),
								new ColumnDefinition(ColumnType.BP_END, Type.LONG),
								new ColumnDefinition(ColumnType.STRAND, Type.STRING),
								new ColumnDefinition(ColumnType.DESCRIPTION, Type.STRING),
								new ColumnDefinition(ColumnType.VALUE, Type.STRING)
						})));
	}

	public GeneParser(FileDefinition fileDefinition) {
		super(fileDefinition);
	}          

	@Override
	public String getName() {
		return "Chipster gene annotation";
	}
	
//	@Override
//	public long getDefaulChunkLength() {
//
//		return 128;
//	}
}