package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class SequenceParser extends TsvParser {

	private static FileDefinition fileDef = new FileDefinition(
			Arrays.asList(
					new ColumnDefinition[] {				
							new ColumnDefinition(ColumnType.CHROMOSOME, Type.LONG),
							new ColumnDefinition(ColumnType.BP_START, Type.LONG),
							new ColumnDefinition(ColumnType.BP_END, Type.LONG),
							new ColumnDefinition(ColumnType.STRAND, Type.STRING),							
							new ColumnDefinition(ColumnType.SEQUENCE, Type.STRING) 
					}));
	
	public SequenceParser() {
		super(fileDef);
	}

	@Override
	public RegionContent[] concise(String chunk) {

		// Return empty table, otherwise TreeNode gets stuck in calling this again
		return new RegionContent[0];
	}

	@Override
	public String getName() {
		return "Chipster sequence file";
	}

//	@Override
//	public long getDefaulChunkLength() {
//		return 1024;
//	}
	
	@Override
	public List<RegionContent> getAll(String chunk, Collection<ColumnType> requestedContents) {

		List<RegionContent> content = super.getAll(chunk, requestedContents);
		System.out.println(content);
		return content;
	}
}