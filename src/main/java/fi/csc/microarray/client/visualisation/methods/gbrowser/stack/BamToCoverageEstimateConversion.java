package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.SortedMap;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.SAMDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ConcisedValueCache;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ConcisedValueCache.Counts;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.CoverageEstimateTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.SamBamUtils;

/**
 * This conversion class uses Picard to read Bam files and creates a coverage estimate by sampling.
 * 
 * @author Aleksi Kallio, Petri Klemelä
 *
 */
public class BamToCoverageEstimateConversion extends SingleThreadAreaRequestHandler {
	
	final static int SAMPLE_SIZE_BP = 1000;

	private SAMDataSource dataSource;
	
	private ConcisedValueCache cache = new ConcisedValueCache();

	public BamToCoverageEstimateConversion(SAMDataSource file, final GBrowser browser) {
	    
		super(null, null);
		
		this.dataSource = file;
	}
		
	@Override
	public void clean() {
		
		SamBamUtils.closeIfPossible(dataSource.getReader());
	}


	@Override
	protected void processAreaRequest(AreaRequest request) {
						
		super.processAreaRequest(request);
		
		if (request.getStatus().poison) {
			return;
		}
			
		processCoverageEstimateRequest(request);									
	}
	
	/**
	 * Return approximation of reads in a given range.
	 * <p>
	 * Works by dividing requested area into several equally-sized regions (samples) and then for each such region asking only a part of its
	 * data (given by SAMPLE_FRACTION) in the beginning of the sample.
	 * <p>
	 * BAM File requests are still quite expensive.
	 * 
	 * @param request
	 * @return
	 */
	public void processCoverageEstimateRequest(AreaRequest request) {

		// How many times file is read
		int step = request.getLength().intValue() / CoverageEstimateTrack.SAMPLING_GRANULARITY;

		int countHits = 0;
		int countReads = 0;
		
		// Divide visible region into subregions and iterate over them
		for (long pos = request.start.bp; pos < request.end.bp; pos += step) {

			BpCoord from = new BpCoord(pos, request.start.chr);
			BpCoord to = new BpCoord(pos + step, request.start.chr);

			// Use cached content, if exists for this region
			SortedMap<BpCoord, Counts> indexedValues = cache.subMap(from, to);
			
			if (!indexedValues.isEmpty()) {
				
				countHits++;
				convertCacheHits(request, step, pos, indexedValues);

			} else {
				countReads++;	
				sampleToGetEstimateRegion(request, from, to);				
			}			
		}
		//System.out.println(countHits + "\t" + countReads);
	}

	private void convertCacheHits(AreaRequest request, int step, long pos, SortedMap<BpCoord, Counts> indexedValues) {

		// Return one result pair for each region covered by one cache hit
		LinkedList<RegionContent> responseList = new LinkedList<RegionContent>(); 
		long cachePos = 0; 
		long startPos = pos;
		int cacheHitsPerRegion = 0;
		for (BpCoord coord : indexedValues.keySet()) {

			// Find end: either next cache hit or end of region
			long endPos;
			SortedMap<BpCoord, Counts> tailMap = indexedValues.tailMap(coord);
			if (tailMap.size() > 1) {
				Iterator<BpCoord> iterator = tailMap.keySet().iterator();
				iterator.next(); // read away this cache hit
				BpCoord next = iterator.next();
				endPos = (startPos + next.bp) / 2;
				cachePos = next.bp;
			} else {
				endPos = (pos + step);
			}
			
			//Region recordRegion = new Region(startPos, endPos, request.start.chr);
			Region recordRegion = new Region(
					cachePos - SAMPLE_SIZE_BP / 2, 
					cachePos + SAMPLE_SIZE_BP / 2, request.start.chr);
			
			LinkedHashMap<ColumnType, Object> values = new LinkedHashMap<ColumnType, Object>();
			values.put(ColumnType.COVERAGE_ESTIMATE_FORWARD, indexedValues.get(coord).forwardCount);
			values.put(ColumnType.COVERAGE_ESTIMATE_REVERSE,  indexedValues.get(coord).reverseCount);
			responseList.add(new RegionContent(recordRegion, values));
			
			cacheHitsPerRegion++;

			// Move to next region
			startPos += endPos;
		}
		
		super.createAreaResult(new AreaResult(request.getStatus(), responseList));
	}

	private void sampleToGetEstimateRegion(AreaRequest request, BpCoord from, BpCoord to) {	
		
		long t = System.currentTimeMillis();
		
		// Fetch new content by taking sample from the middle of this area
		long stepMiddlepoint = (from.bp + to.bp) / 2;
		long start = stepMiddlepoint - SAMPLE_SIZE_BP / 2;
		long end = stepMiddlepoint + SAMPLE_SIZE_BP / 2;
		CloseableIterator<SAMRecord> iterator = dataSource.getReader().query(dataSource.getChromosomeNameUnnormaliser().unnormalise(from.chr), (int) start, (int) end, false);

		// Count reads in this sample area
		int countForward = 0;
		int countReverse = 0;

		int countRejected1 = 0;
		int countRejected2 = 0;
		
		
		for (Iterator<SAMRecord> i = iterator; i.hasNext();) {			
			
			if (super.hasNewRequest()) {
				//Stop iteration, but send results
				break;
			}
			
			SAMRecord record = i.next();

			// Accept only records that start in this area (very rough approximation for spliced reads)
			if (record.getAlignmentStart() >= start) { 
				if (record.getAlignmentEnd() <= end) {


					if (record.getReadNegativeStrandFlag()) {
						countReverse++;
					} else {
						countForward++;
					}
				} else {
					countRejected2++;
					 
					//Stop iteration, because Picard may return many reads after the requested reqion  
					break;
				}
			}else {
				countRejected1++;
			}
		}

		//System.out.println("Forward: " + countForward + "\tReverse: " + countReverse + "\tRejected1: " + countRejected1 + "\tRejected2: " + countRejected2 + "\t" + (System.currentTimeMillis() - t) + " ms");
		
		
		cache.store(new BpCoord(stepMiddlepoint, from.chr), countForward, countReverse);
		
		// We are done
		iterator.close();

		// Send result
		LinkedList<RegionContent> content = new LinkedList<RegionContent>();
		
		LinkedHashMap<ColumnType, Object> values = new LinkedHashMap<ColumnType, Object>();
		values.put(ColumnType.COVERAGE_ESTIMATE_FORWARD, countForward);
		values.put(ColumnType.COVERAGE_ESTIMATE_REVERSE, countReverse);
		
		content.add(new RegionContent(new Region(start, end, from.chr), values));
		
		super.createAreaResult(new AreaResult(request.getStatus(), content));						
	}
	
	public String toString() {
		return this.getClass().getName() + " - " + dataSource;
	}
}
