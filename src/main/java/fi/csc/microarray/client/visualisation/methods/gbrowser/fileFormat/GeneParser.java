package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class GeneParser extends TsvParser{

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
	public RegionContent[] concise(BpCoordRegion nodeBpRegion) {

		long totalF = 0;
		long totalR = 0;
		
		String[] first = getFirstRow();
		String[] last = getLastRow();
		
		BpCoord start = new BpCoord((Long)get(first, ColumnType.BP_START), (Chromosome)get(first, ColumnType.CHROMOSOME));
		BpCoord end = new BpCoord((Long)get(last, ColumnType.BP_START), (Chromosome)get(last, ColumnType.CHROMOSOME));		
		BpCoordRegion reg = new BpCoordRegion(start, end);		
		long length =  reg.getLength();
		
		for (RegionContent rc : 
			getAll(Arrays.asList(new ColumnType[] { ColumnType.STRAND }))) {

			if((Strand)rc.values.get(ColumnType.STRAND) == Strand.FORWARD) {
				
				totalF += length;
				
			} else {
				totalR += length;
			}			
		}
		
		RegionContent[] result = new RegionContent[] {
			new RegionContent(nodeBpRegion, totalF / (float)nodeBpRegion.getLength()),
			new RegionContent(nodeBpRegion, totalR / (float)nodeBpRegion.getLength())
		};		
		
		result[0].values.put(ColumnType.STRAND, Strand.FORWARD);
		result[1].values.put(ColumnType.STRAND, Strand.REVERSED);
		
		return result;
	}

	@Override
	public FileParser clone() {
		FileParser clone = new GeneParser();
		
		clone.chunk = this.chunk;
		
		return clone;
	}

	@Override
	public String getName() {
		return "Chipster gene annotation";
	}
}