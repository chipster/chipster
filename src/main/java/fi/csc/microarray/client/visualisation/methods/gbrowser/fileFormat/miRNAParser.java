package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;

public class miRNAParser extends GeneParser{

	public miRNAParser() {
		super(new FileDefinition(
				Arrays.asList(
						new ColumnDefinition[] {

								new ColumnDefinition(ColumnType.CHROMOSOME, Type.LONG),
								new ColumnDefinition(ColumnType.BP_START, Type.LONG),
								new ColumnDefinition(ColumnType.BP_END, Type.LONG),
								new ColumnDefinition(ColumnType.STRAND, Type.STRING),
								new ColumnDefinition(ColumnType.DESCRIPTION, Type.STRING)

						})));
	}
	
	@Override
	public FileParser clone() {
		FileParser clone = new miRNAParser();	
		clone.chunk = this.chunk;		
		return clone;
	}
}