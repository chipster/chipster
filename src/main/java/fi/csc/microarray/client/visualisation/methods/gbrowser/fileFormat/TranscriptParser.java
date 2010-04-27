package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

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
								new ColumnDefinition(ColumnType.PARENT_PART, Type.LONG) // should be Type.LONG, but there are some quotes in current annotation file
		})));
	}

	@Override
	public FileParser clone() {
		FileParser clone = new TranscriptParser();
		clone.chunk = (chunk == null ? null : new String(this.chunk));
		return clone;
	}

	@Override
	public String getName() {
		return "Chipster transcript file";
	}
	
	@Override
	public Object get(String[] cols, ColumnType col) {
		
		return super.get(cols, col);
	}

	@Override
	public List<RegionContent> getAll(Collection<ColumnType> requestedContents) {

		return super.getAll(requestedContents); 
	}
}