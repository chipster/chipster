package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class HeaderTsvParser extends TsvParser {

	public HeaderTsvParser() {
		super(new FileDefinition(Arrays.asList(
				new ColumnDefinition[] { 
						new ColumnDefinition(ColumnType.CHROMOSOME, Type.STRING), 
						new ColumnDefinition(ColumnType.BP_START, Type.LONG), 
						new ColumnDefinition(ColumnType.BP_END, Type.LONG), 
						new ColumnDefinition(ColumnType.SKIP, Type.STRING), 
						new ColumnDefinition(ColumnType.SKIP, Type.STRING), 
						new ColumnDefinition(ColumnType.SKIP, Type.STRING), 
						new ColumnDefinition(ColumnType.SKIP, Type.STRING), 
				})));
	}

	public HeaderTsvParser(FileDefinition fileDefinition) {
		super(fileDefinition);
	}

	@Override
	public String[] getFirstRow(String chunk) {
		return super.getFirstRow(filterChunk(chunk));
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
	public RegionContent[] concise(String chunk) {
		return new RegionContent[] {};
	}

	@Override
	public Object get(String[] cols, ColumnType col) {

		Object obj = super.get(cols, col);

		if (col == ColumnType.CHROMOSOME) {
			return new Chromosome(((Chromosome) obj).toString().replace(".fa", ""));
		}
		
		return obj;
	}
	
	@Override
	public List<RegionContent> getAll(String chunk, Collection<ColumnType> requestedContents) {
		return super.getAll(filterChunk(chunk), requestedContents);
	}

	@Override
	public BpCoordRegion getBpRegion(String chunk) {
		return super.getBpRegion(filterChunk(chunk));
	};
	
	@Override
	public String[] getLastRow(String chunk) {
		return super.getLastRow(filterChunk(chunk));
	};
	
	
	private String filterChunk(String chunk) {
		String newChunk = "";
		for (String row : chunk.split("\n")) {
			if (!row.startsWith("chr\tstart")) {
				newChunk += (row + "\n");
			}
		}
		return newChunk;
	}
}