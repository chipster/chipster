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
import net.sf.picard.reference.ChipsterIndexedFastaSequenceFile;
import net.sf.picard.reference.ReferenceSequence;
import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.IndexedFastaDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.LineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;


public class IndexedFastaFileFetcherThread extends Thread {

	private BlockingQueue<BpCoordFileRequest> fileRequestQueue;
	private ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue;

	private IndexedFastaDataSource dataSource;

	private IndexedFastaHandlerThread areaRequestThread;

	private boolean poison = false;

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


		// Process only new part of the requested area
		// FIXME We rely on other layers to cache previous results, which is not very clean.
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

		ChunkDataSource data = new ChunkDataSource(dataSource.getUrl(), null, null);
		LineDataSource index = new LineDataSource(dataSource.getIndex(), null);

		ChipsterIndexedFastaSequenceFile picard = new ChipsterIndexedFastaSequenceFile(data, index);

		List<RegionContent> responseList = new LinkedList<RegionContent>();

		try {
			ReferenceSequence seq = picard.getSubsequenceAt(request.start.chr.toNormalisedString(), request.start.bp, request.end.bp);


			Region recordRegion = new Region(request.start.bp, request.end.bp, request.start.chr);

			LinkedHashMap<ColumnType, Object> values = new LinkedHashMap<ColumnType, Object>();

			RegionContent regCont = new RegionContent(recordRegion, values);

			values.put(ColumnType.SEQUENCE, new String(seq.getBases()));

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
