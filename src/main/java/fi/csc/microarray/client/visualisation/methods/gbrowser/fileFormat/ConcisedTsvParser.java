package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import com.sun.xml.internal.messaging.saaj.packaging.mime.internet.ContentType;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public abstract class ConcisedTsvParser extends TsvParser{

	public ConcisedTsvParser(FileDefinition fileDef) {
		super(fileDef);
	}
	
	@Override
	public RegionContent[] concise(String chunk) {
		
		Long start = (Long)get(getFirstRow(chunk), ColumnType.BP_START);
		Long end = (Long)get(getLastRow(chunk), ColumnType.BP_START);
		Chromosome startChr = (Chromosome)get(getFirstRow(chunk), ColumnType.CHROMOSOME);
		Chromosome endChr = (Chromosome)get(getLastRow(chunk), ColumnType.CHROMOSOME);
				
		Long readLength =  (new BpCoordRegion(start, end, startChr)).getLength();
		
		final long binSize = 1000;
		
		long firstBin = (long)Math.floor(start / binSize) * binSize;
		long lastBin = (long)Math.floor(end / binSize) * binSize;
		
		int binCount = (int) ((lastBin - firstBin) / binSize + 1);
		

		if(readLength != null && startChr.equals(endChr) && binCount > 0) {
			
			long[] fBins = new long[binCount];
			long[] rBins = new long[binCount];

			for (RegionContent rc : 
				getAll(chunk, Arrays.asList(new ColumnType[] { ColumnType.STRAND, ColumnType.BP_START }))) {
				
				long readStart = (Long)rc.values.get(ColumnType.BP_START);
				int bin = (int) ((readStart - firstBin) / binSize);
				
				if((Strand)rc.values.get(ColumnType.STRAND) == Strand.FORWARD) {

					fBins[bin] += readLength;

				} else {
					rBins[bin] += readLength;
				}
			}
			
			List<RegionContent> results = new LinkedList<RegionContent>();
			
			for (int i = 0; i < binCount; i ++) {
				
				long binStart = firstBin + i * binSize;
				long binEnd = binStart + binSize - 1;
				
				RegionContent fRc = new RegionContent(new BpCoordRegion(binStart, binEnd, startChr), fBins[i] / (float)binSize);
				RegionContent rRc = new RegionContent(new BpCoordRegion(binStart, binEnd, startChr), rBins[i] / (float)binSize);
				
				fRc.values.put(ColumnType.STRAND, Strand.FORWARD);
				rRc.values.put(ColumnType.STRAND, Strand.REVERSED);
				
				//System.out.println(fRc.region + ", " + fRc.values.get(ColumnType.VALUE));
				
				results.add(fRc);
				results.add(rRc);
			}

			return results.toArray(new RegionContent[results.size()]);
		} else {
			
			//FIXME Length of region can't be calculated, because it contains two or more chromosomes, do something wise
			return new  RegionContent[] {};
		}
	}
}
