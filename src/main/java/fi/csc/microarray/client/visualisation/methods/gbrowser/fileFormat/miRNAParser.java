package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;

public class miRNAParser extends GeneParser{

	public miRNAParser() {
		super(new FileDefinition(
				Arrays.asList(
						new ColumnDefinition[] {

								new ColumnDefinition(ColumnType.CHROMOSOME, Type.LONG, 2),
								new ColumnDefinition(ColumnType.BP_START, Type.LONG, 16),
								new ColumnDefinition(ColumnType.BP_END, Type.LONG, 16),
								new ColumnDefinition(ColumnType.STRAND, Type.STRING, 2),
								new ColumnDefinition(ColumnType.DESCRIPTION, Type.STRING, 32),					
								new ColumnDefinition(ColumnType.SKIP, Type.NEWLINE, 1)

						})));
	}
	
	@Override
	public FileParser clone() {
		FileParser clone = new miRNAParser();
		
		clone.chunk = this.chunk;
		
		return clone;
	}
}