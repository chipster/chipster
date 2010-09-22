package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class SNPParser extends TsvParser{

	public SNPParser() {
		super(new FileDefinition(
				Arrays.asList(
						new ColumnDefinition[] {
								new ColumnDefinition(ColumnType.ID, Type.STRING),
								new ColumnDefinition(ColumnType.CHROMOSOME, Type.LONG),
								new ColumnDefinition(ColumnType.POSITION, Type.LONG), //position on chromosome
								new ColumnDefinition(ColumnType.ALLELE, Type.STRING),
								new ColumnDefinition(ColumnType.STRAND, Type.STRING),
								new ColumnDefinition(ColumnType.CONSEQUENCE_TO_TRANSCRIPT, Type.STRING)
						})));
	}
	
	@Override
	public RegionContent[] concise(Chunk chunk) {
		return null;
	}

	@Override
	public String getName() {
		return "SNP Parser";
	}
	
	@Override
	public BpCoordRegion getBpRegion(Chunk chunk) {
		String[] firstRow = getFirstRow(chunk);
	    String[] lastRow = getLastRow(chunk);
		Long start = (Long)get(firstRow, ColumnType.POSITION);
		Long end;
        end = (Long)get(lastRow, ColumnType.POSITION);
        
		Chromosome startChr = (Chromosome)get(firstRow, ColumnType.CHROMOSOME);
		Chromosome endChr = (Chromosome)get(lastRow, ColumnType.CHROMOSOME);
		return new BpCoordRegion(start, startChr, end, endChr);
	}
	
	@Override
	public List<RegionContent> getAll(Chunk chunk,
			Collection<ColumnType> requestedContents) {
		List<RegionContent> rows = new LinkedList<RegionContent>();
		
		for (String row : chunk.getContent().split("\n")) {
			
			Map<ColumnType, Object> values = new HashMap<ColumnType, Object>();
			
			String[] cols = row.split("\t");
			for (ColumnType requestedContent : requestedContents) {
				values.put(requestedContent, this.get(cols, requestedContent));					
			}
			
			Long start = (Long)get(cols, ColumnType.POSITION);
			Long end = (Long)get(cols, ColumnType.POSITION);
			Chromosome chr = (Chromosome)get(cols, ColumnType.CHROMOSOME);
	
			rows.add(new RegionContent(new BpCoordRegion(start, end, chr), values));

		}
		return rows;
	}

	@Override
	public long getHeaderLength(File file) {
//		try {
//			FileReader fileReader = new FileReader(file);
//			BufferedReader reader = new BufferedReader(fileReader);
//			try {
//				System.out.println(reader.readLine().length());
//				return reader.readLine().length();
//			} catch (IOException e) {
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		} catch (FileNotFoundException e) {
//			e.printStackTrace();
//		}
		return 0;
	}
}
