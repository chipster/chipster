package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Subclass of {@link TsvParser} implementing the concised method for our own annotation file format.
 * 
 * @author Petri Klemelä
 *
 */
public abstract class AnnotationTsvParser extends TsvParser{

	public AnnotationTsvParser(FileDefinition fileDef) {
		super(fileDef);
	}
	
	@Override
	public RegionContent[] concise(Chunk chunk) {

		long totalF = 0;
		long totalR = 0;
		
		Long start = (Long)get(getFirstRow(chunk), ColumnType.BP_START);
		Long end = (Long)get(getLastRow(chunk), ColumnType.BP_START);
		Chromosome chr = (Chromosome)get(getFirstRow(chunk), ColumnType.CHROMOSOME);

		Long readLength =  (new BpCoordRegion(start, end, chr)).getLength();

		if(readLength != null) {

			for (RegionContent rc : 
				getAll(chunk, Arrays.asList(new ColumnType[] { ColumnType.STRAND }))) {

				if((Strand)rc.values.get(ColumnType.STRAND) == Strand.FORWARD) {

					totalF += readLength;

				} else {
					totalR += readLength;
				}
			}
			
			Long regionLength = getBpRegion(chunk).getLength();

			if (regionLength != null) {
				RegionContent[] result = new RegionContent[] {
						new RegionContent(getBpRegion(chunk), totalF / (float)regionLength),
						new RegionContent(getBpRegion(chunk), totalR / (float)regionLength)
				};		

				result[0].values.put(ColumnType.STRAND, Strand.FORWARD);
				result[1].values.put(ColumnType.STRAND, Strand.REVERSED);

				return result;
			}
		} 			
		//FIXME Length of region can't be calculated, because it contains two or more chromosomes, do something wise
		return new  RegionContent[] {};
		
	}
}