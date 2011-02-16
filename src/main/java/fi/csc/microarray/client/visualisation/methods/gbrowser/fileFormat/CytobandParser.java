package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Parser for cytoband description files using tab-separated values.
 *
 */
public class CytobandParser extends TsvParser {

	public static final String LAST_ROW_OF_CHROMOSOME = "lastRowOfChromosome";
	public static final String TRUE = "true";

	public CytobandParser() {
		super(new FileDefinition(
				Arrays.asList(
						new ColumnDefinition[] {
								new ColumnDefinition(ColumnType.CHROMOSOME, Type.STRING),
								new ColumnDefinition(ColumnType.BP_START, Type.LONG),
								new ColumnDefinition(ColumnType.BP_END, Type.LONG),
								new ColumnDefinition(ColumnType.ID, Type.STRING),
								new ColumnDefinition(ColumnType.VALUE, Type.STRING)
						})));
	}

	@Override
	public RegionContent[] concise(Chunk chunk) {
		// return empty table, otherwise TreeNode gets stuck in calling this again
		return new RegionContent[0];
	}

	@Override
	public String getName() {
		return "Cytobands";
	}
	
	/**
	 * Overridden to notice the chromosome changes. This information is needed to find out the last
	 * part of each chromosome to set the scroll limits. If the chromosome change happens to 
	 * be between chunks, it won't be noticed and that chromosome doesn't get limited.
	 * 
	 * @see fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.TsvParser#getAll(java.lang.String, java.util.Collection)
	 */
	@Override
	public List<RegionContent> getAll(Chunk chunk, Collection<ColumnType> requestedContents) {

		List<RegionContent> rows = new LinkedList<RegionContent>();
		

		Chromosome lastChr = null;
		
		for (String row : chunk.getContent().split("\n")) {
			
			LinkedHashMap<ColumnType, Object> values = new LinkedHashMap<ColumnType, Object>();
			
			String[] cols = row.split("\t");
			
			for (ColumnType requestedContent : requestedContents) {
						
				values.put(requestedContent, this.get(cols, requestedContent));					
			}
			
			Long start = (Long)get(cols, ColumnType.BP_START);
			Long end = (Long)get(cols, ColumnType.BP_END);
			Chromosome chr = (Chromosome)get(cols, ColumnType.CHROMOSOME);
			
			if (lastChr != null && !lastChr.equals(chr)) {
				
				Map<String, String> metaMap = new HashMap<String, String>();
				metaMap.put(LAST_ROW_OF_CHROMOSOME, TRUE);
				
				rows.get(rows.size() - 1).values.put(ColumnType.METADATA, metaMap);
			}
				
			rows.add(new RegionContent(new BpCoordRegion(start, end, chr), values));
			
			lastChr = chr;
		}
		
		return rows;
	}
}