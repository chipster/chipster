package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.SortedMap;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.LinkedBlockingQueue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.SAMDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ConcisedValueCache;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ConcisedValueCache.Counts;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.IntensityTrack;

/**
 * The processing layer thread for SAM compatible files.
 * 
 * @author Aleksi Kallio
 *
 */
public class SAMHandlerThread extends AreaRequestHandler {
    
	private ConcisedValueCache cache = new ConcisedValueCache();

	private SAMDataSource samData;
	private SAMFileFetcherThread fileFetcher;
	private BlockingQueue<BpCoordFileRequest> fileRequestQueue = new LinkedBlockingQueue<BpCoordFileRequest>();
	private ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue = new ConcurrentLinkedQueue<ParsedFileResult>();

    public SAMHandlerThread(DataSource file, Queue<AreaRequest> areaRequestQueue,
            AreaResultListener areaResultListener) {
        
        super(areaRequestQueue, areaResultListener);
        samData = (SAMDataSource) file;
    }

	@Override
	public synchronized void run() {

		// Start file processing layer thread
		fileFetcher = new SAMFileFetcherThread(fileRequestQueue, fileResultQueue, this, samData);
		fileFetcher.start();
		
		// Start this thread
		super.run();
	}

	protected boolean checkOtherQueues() {
		ParsedFileResult fileResult = null;
		if ((fileResult = fileResultQueue.poll()) != null) {
			processFileResult(fileResult);
		}
		return fileResult != null;
	}

    private void processFileResult(ParsedFileResult fileResult) {

    	if (fileResult.getStatus().concise) {

    		LinkedList<RegionContent> responseList = new LinkedList<RegionContent>();
    		for (RegionContent content : fileResult.getContents()) {
    			int countForward = (Integer) content.values.get(ColumnType.VALUE_FORWARD); 
    			int countReverse = (Integer) content.values.get(ColumnType.VALUE_REVERSE);
    			cache.store(new BpCoord(content.region.getMid(), content.region.start.chr), countForward, countReverse);
    			addConcisedRegionContents(fileResult.areaRequest, responseList, content.region.start.bp, content.region.end.bp, countForward, countReverse);
    		}
    		
    		createAreaResult(new AreaResult(fileResult.getStatus(), responseList));
    		
    	} else {
    		createAreaResult(new AreaResult(fileResult.getStatus(), fileResult.getContents()));
    	}
	}

	/**
     * Handles normal and concised area requests by using SAMFile.
     */
    @Override
    protected void processAreaRequest(AreaRequest areaRequest) {
		if (areaRequest.status.concise) {
			processConcisedAreaRequest(areaRequest);

		} else {
			fileRequestQueue.add(new BpCoordFileRequest(areaRequest, areaRequest.start, areaRequest.end, areaRequest.status));
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
	 */
	public void processConcisedAreaRequest(AreaRequest request) {

		// How many times file is read
		int step = request.getLength().intValue() / IntensityTrack.SAMPLING_GRANULARITY;

		// Divide visible region into subregions and iterate over them
		for (long pos = request.start.bp; pos < request.end.bp; pos += step) {

			BpCoord from = new BpCoord(pos, request.start.chr);
			BpCoord to = new BpCoord(pos + step, request.start.chr);

			// Use cached content, if exists for this region
			SortedMap<BpCoord, Counts> indexedValues = cache.subMap(from, to);
			if (!indexedValues.isEmpty()) {
				convertCacheHitsToConcisedRegions(request, step, pos, indexedValues);

			} else {
				
				request.status.maybeClearQueue(fileRequestQueue);
				
				fileRequestQueue.add(new BpCoordFileRequest(request, from, to, request.status));
			}

		}
	}

	private void convertCacheHitsToConcisedRegions(AreaRequest request, int step, long pos, SortedMap<BpCoord, Counts> indexedValues) {

		// Return one result pair for each region covered by one cache hit
		LinkedList<RegionContent> responseList = new LinkedList<RegionContent>(); 
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

			} else {
				endPos = (pos + step);
			}

			addConcisedRegionContents(request, responseList, startPos, endPos, indexedValues.get(coord).forwardCount, indexedValues.get(coord).reverseCount);
			cacheHitsPerRegion++;

			// Move to next region
			startPos += endPos;
		}
		
		createAreaResult(new AreaResult(request.status, responseList));
	}

	private void addConcisedRegionContents(AreaRequest request, List<RegionContent> responseList, long startPos, long endPos, int countForward, int countReverse) {
		// Create two approximated response objects: one for each strand
		BpCoordRegion recordRegion = new BpCoordRegion(startPos, endPos, request.start.chr);

		// Forward
		LinkedHashMap<ColumnType, Object> values = new LinkedHashMap<ColumnType, Object>();
		values.put(ColumnType.VALUE, (float) countForward);
		values.put(ColumnType.STRAND, Strand.FORWARD);
		responseList.add(new RegionContent(recordRegion, values));

		// Reverse
		values = new LinkedHashMap<ColumnType, Object>();
		values.put(ColumnType.VALUE, (float) countReverse);
		values.put(ColumnType.STRAND, Strand.REVERSED);
		responseList.add(new RegionContent(recordRegion, values));
	}

}
