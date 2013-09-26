package fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex;

import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.SortedMap;

import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex.ConcisedValueCache.Counts;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.CoverageEstimateTrack;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;

/**
 * This conversion class uses Picard to read Bam files and creates a coverage estimate by sampling.
 * 
 * @author Aleksi Kallio, Petri Klemelä
 *
 */
public class BamToCoverageEstimateConversion extends DataThread {
	
	public final static int SAMPLE_SIZE_BP = 1000;

	private BamDataSource dataSource;
	
	private ConcisedValueCache cache = new ConcisedValueCache();

	public BamToCoverageEstimateConversion(BamDataSource file, final GBrowser browser) {
	    
		super(browser, file);
		
		this.dataSource = file;
	}
		
	@Override
	public void clean() {
		
		dataSource.close();
	}


	@Override
	protected void processDataRequest(DataRequest request) throws GBrowserException {					
		
		if (request.getRequestedContents().contains(DataType.CANCEL)) {
			return;
		}
		
		processCoverageEstimateRequest(request);
		
		if (!hasNewRequest()) {

			createBetterEstimate(request, 4);
		}
		
		if (!hasNewRequest()) {

			createBetterEstimate(request, 16);
		}
	}
	
	private void createBetterEstimate(DataRequest request, int partCount) throws GBrowserException {
		
		long step = getDataRegion().getLength() / partCount;
		
		if (step < SAMPLE_SIZE_BP) {
			//no need for better estimate
			return;
		}
		
		for (long i = getDataRegion().start.bp; i < getDataRegion().end.bp - step; i += step) {
			
			Region region = new Region(i, Math.min(i + step, getDataRegion().end.bp), getDataRegion().start.chr);
			DataRequest requestPart;

			requestPart = new DataRequest(region, request.getRequestedContents(), new DataStatus(request.getStatus()));

			if (this.hasNewRequest()) {
				return;
			} else {
				processCoverageEstimateRequest(requestPart);
			}			
		}		
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
	 * @throws GBrowserException 
	 */
	public void processCoverageEstimateRequest(DataRequest request) throws GBrowserException {

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
		//System.out.println("Cache: " + countHits + "\tFile: " + countReads);
	}

	private void convertCacheHits(DataRequest request, int step, long pos, SortedMap<BpCoord, Counts> indexedValues) {

		// Return one result pair for each region covered by one cache hit
		LinkedList<Feature> responseList = new LinkedList<Feature>(); 
		long cachePos = 0; 

		for (BpCoord coord : indexedValues.keySet()) {

			cachePos = coord.bp;
			
			Region recordRegion = new Region(
					cachePos - SAMPLE_SIZE_BP / 2, 
					cachePos + SAMPLE_SIZE_BP / 2, request.start.chr);
			
			LinkedHashMap<DataType, Object> values = new LinkedHashMap<DataType, Object>();
			values.put(DataType.COVERAGE_ESTIMATE_FORWARD, indexedValues.get(coord).forwardCount);
			values.put(DataType.COVERAGE_ESTIMATE_REVERSE,  indexedValues.get(coord).reverseCount);
			responseList.add(new Feature(recordRegion, values));
		}
		
		super.createDataResult(new DataResult(request.getStatus(), responseList));
	}

	private void sampleToGetEstimateRegion(DataRequest request, BpCoord from, BpCoord to) throws GBrowserException {	
		
		long t = System.currentTimeMillis();
		
		// Fetch new content by taking sample from the middle of this area
		long stepMiddlepoint = (from.bp + to.bp) / 2;
		long start = stepMiddlepoint - SAMPLE_SIZE_BP / 2;
		long end = stepMiddlepoint + SAMPLE_SIZE_BP / 2;
		
		CloseableIterator<SAMRecord> iterator = null;
		
		try {
			iterator = dataSource.query(from.chr, (int) start, (int) end);
		
		} catch (RuntimeException e) {
			throw new GBrowserException("Error in data query. Ususally this happens when a wrong index file is selected.", e);
		} catch (OutOfMemoryError e) {
			//Not a real outOfMemory, picard just tried to create an arbitrarily large buffer 
			throw new GBrowserException("Error in data query. Ususally this happens when a wrong index file is selected. " + e.getMessage());
		}

		// Count reads in this sample area
		int countForward = 0;
		int countReverse = 0;

		int countRejected1 = 0;
		int countRejected2 = 0;
		
		boolean interrupted = false;
		
		for (Iterator<SAMRecord> i = iterator; i.hasNext();) {			
			
//			if (super.hasNewRequest()) {
//				//Stop iteration, but send results
//				interrupted = true;
//				break;
//			}
			
			SAMRecord record = i.next();

			// Accept only records that start in this area (very rough approximation for spliced reads)
			if (record.getAlignmentStart() >= start) { 
				if (record.getAlignmentEnd() <= end) {


					if (record.getReadNegativeStrandFlag()) {
						countReverse += record.getReadLength();
					} else {
						countForward += record.getReadLength();
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
		
		// We are done
		iterator.close();
		
		if (!interrupted) {
			cache.store(new BpCoord(stepMiddlepoint, from.chr), countForward, countReverse);		

			// Send result
			LinkedList<Feature> content = new LinkedList<Feature>();

			LinkedHashMap<DataType, Object> values = new LinkedHashMap<DataType, Object>();
			values.put(DataType.COVERAGE_ESTIMATE_FORWARD, countForward);
			values.put(DataType.COVERAGE_ESTIMATE_REVERSE, countReverse);

			content.add(new Feature(new Region(start, end, from.chr), values));

			super.createDataResult(new DataResult(request.getStatus(), content));
		}
	}
	
	public String toString() {
		return this.getClass().getName() + " - " + dataSource;
	}
}
