package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;


import java.util.Arrays;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class GeneIndexParser extends TsvParser{

	public GeneIndexParser() {
		super(new FileDefinition(
				Arrays.asList(
						new ColumnDefinition[] {
								new ColumnDefinition(ColumnType.CHROMOSOME, Type.LONG),
								new ColumnDefinition(ColumnType.BP_START, Type.LONG),
								new ColumnDefinition(ColumnType.BP_END, Type.LONG),
								new ColumnDefinition(ColumnType.SKIP, Type.STRING),
								new ColumnDefinition(ColumnType.DESCRIPTION, Type.STRING),
								new ColumnDefinition(ColumnType.SKIP, Type.STRING)
						})));
	}

	@Override
	public RegionContent[] concise(String chunk) {
		return null;
	}

	@Override
	public String getName() {
		return "Gene indexes parser";
	}

}
