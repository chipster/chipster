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

public class PeakParser extends TsvParser {

	public PeakParser() {
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

	public PeakParser(FileDefinition fileDefinition) {
		super(fileDefinition);
	}

	@Override
	public String[] getFirstRow(String chunk) {
		String row = chunk.substring(0, chunk.indexOf("\n"));
		if (row.startsWith("chr")) {
			row = chunk.substring(chunk.indexOf("\n") + 1, chunk.indexOf("\n", chunk.indexOf("\n") + 1));
		}
		return row.split("\t");
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

		List<RegionContent> rows = new LinkedList<RegionContent>();
	
		for (String row : chunk.split("\n")) {
	
			if (row.startsWith("chr")) {
				continue;
			}
			
			Map<ColumnType, Object> values = new HashMap<ColumnType, Object>();
			
			String[] cols = row.split("\t");
			
			for (ColumnType requestedContent : requestedContents) {
						
				values.put(requestedContent, this.get(cols, requestedContent));					
			}
			
			Long start = (Long)get(cols, ColumnType.BP_START);
			Long end = (Long)get(cols, ColumnType.BP_END);
			Chromosome chr = (Chromosome)get(cols, ColumnType.CHROMOSOME);
	
			rows.add(new RegionContent(new BpCoordRegion(start, end, chr), values));

		}
		
		return rows;
	}
}