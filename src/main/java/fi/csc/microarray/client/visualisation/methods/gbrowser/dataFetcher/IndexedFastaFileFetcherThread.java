package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URISyntaxException;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;

import net.sf.picard.PicardException;
import net.sf.picard.reference.ReferenceSequence;
import fi.csc.microarray.client.visualisation.methods.gbrowser.IndexedFastaDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;


public class IndexedFastaFileFetcherThread extends Thread {

	private static final Long BUFFER_EXTRA = 1000l;
	private BlockingQueue<BpCoordFileRequest> fileRequestQueue;
	private ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue;

	private IndexedFastaDataSource dataSource;

	private IndexedFastaHandlerThread areaRequestThread;

	private boolean poison = false;
	private SequenceBuffer buffer;

	public IndexedFastaFileFetcherThread(BlockingQueue<BpCoordFileRequest> fileRequestQueue, ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue, IndexedFastaHandlerThread areaRequestThread, IndexedFastaDataSource dataSource) {

		this.fileRequestQueue = fileRequestQueue;
		this.fileResultQueue = fileResultQueue;
		this.areaRequestThread = areaRequestThread;
		this.dataSource = dataSource;
		this.setDaemon(true);
	}

	public void run() {

		while (!poison) {
			try {

				for (BpCoordFileRequest fileRequest : fileRequestQueue) {
					if (fileRequest.getStatus().poison) {
						poison = true;
						return;
					}
				}

				processFileRequest(fileRequestQueue.take());

			} catch (IOException e) {
				e.printStackTrace(); // FIXME fix exception handling
			} catch (InterruptedException e) {
				e.printStackTrace(); // FIXME fix exception handling
			} catch (URISyntaxException e) {
				e.printStackTrace();
			}
		}
	}

	private void processFileRequest(BpCoordFileRequest fileRequest) throws IOException, URISyntaxException {

		if (fileRequest.getStatus().poison) {
			poison = true;
			return;
		}


		AreaRequest request = fileRequest.areaRequest;

		fetchSequence(fileRequest);
	}

	/**
	 * Find reads in a given range.
	 * 
	 * <p>
	 * TODO add cigar to the list of returned values
	 * <p>
	 * TODO add pair information to the list of returned values
	 * 
	 * @param request
	 * @return
	 * @throws URISyntaxException 
	 * @throws FileNotFoundException 
	 */
	public void fetchSequence(BpCoordFileRequest fileRequest) throws FileNotFoundException, URISyntaxException {

		AreaRequest request = fileRequest.areaRequest;

		List<RegionContent> responseList = new LinkedList<RegionContent>();

		try {
			
			String chr = dataSource.getChromosomeNameUnnormaliser().unnormalise(request.start.chr);

			//A few thousand byte buffer eliminates the need for downloading data on every frame. 
			//Buffer has to be refreshed relatively rarely, maybe in 1 second intervals during heavy scrolling
			if (buffer == null  || !buffer.contains(request.start.bp, request.end.bp, request.start.chr)) {

				long bufferStart = Math.max(request.start.bp - BUFFER_EXTRA, 0);
				long bufferEnd = request.end.bp + BUFFER_EXTRA;				
				
				ReferenceSequence seq = dataSource.getPicard().getSubsequenceAt(chr, bufferStart, bufferEnd);
				
				buffer = new SequenceBuffer(new String(seq.getBases()), bufferStart, request.start.chr);
			}

			String seqString = buffer.get(request.start.bp, request.end.bp, request.start.chr);

			Region recordRegion = new Region(request.start.bp, request.end.bp, request.start.chr);

			LinkedHashMap<ColumnType, Object> values = new LinkedHashMap<ColumnType, Object>();
			values.put(ColumnType.SEQUENCE, seqString);

			RegionContent regCont = new RegionContent(recordRegion, values);


			/*
			 * NOTE! RegionContents created from the same read area has to be equal in methods equals, hash and compareTo. Primary types
			 * should be ok, but objects (including tables) has to be handled in those methods separately. Otherwise tracks keep adding
			 * the same reads to their read sets again and again.
			 */
			responseList.add(regCont);

		} catch (PicardException e) {
			e.printStackTrace(); //Catch "Query asks for data past end of contig" to prevent this thread from ending
		}

		// Send result
		ParsedFileResult result = new ParsedFileResult(responseList, fileRequest, fileRequest.areaRequest, fileRequest.getStatus());
		fileResultQueue.add(result);
		areaRequestThread.notifyAreaRequestHandler();
	}

	public String toString() {
		return this.getClass().getName() + " - " + dataSource;
	}
}
