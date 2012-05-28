package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.Queue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.LinkedBlockingQueue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.LineDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneRequest;
import fi.csc.microarray.gbrowser.index.GeneResult;

/**
 * Processing layer of gene and transcript annotation data. There is no need for processing with this data type, so just passing 
 * original data through.
 * 
 * @author Petri Klemelä
 *
 */
public class GtfHandlerThread extends AreaRequestHandler {

	private LineDataSource dataSource;

	private GtfFileFetcherThread fileFetcher;
	private BlockingQueue<BpCoordFileRequest> fileRequestQueue = new LinkedBlockingQueue<BpCoordFileRequest>();
	private ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue = new ConcurrentLinkedQueue<ParsedFileResult>();

	public GtfHandlerThread(DataSource file, Queue<AreaRequest> areaRequestQueue,
			AreaResultListener areaResultListener) {

		super(areaRequestQueue, areaResultListener);
		dataSource = (LineDataSource) file;
	}

	@Override
	public synchronized void run() {

		// Start file processing layer thread
		fileFetcher = new GtfFileFetcherThread(fileRequestQueue, fileResultQueue, this,
				dataSource);

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

		if (fileResult.getFileRequest().areaRequest instanceof GeneRequest) {
			GeneRequest geneRequest = (GeneRequest)fileResult.getFileRequest().areaRequest;
			createAreaResult(new GeneResult(fileResult.getStatus(), fileResult.getContents(), geneRequest.getSearchString()));
		} else {
			createAreaResult(new AreaResult(fileResult.getStatus(), fileResult.getContents()));
		}
	}

	@Override
	protected void processAreaRequest(AreaRequest areaRequest) {

		super.processAreaRequest(areaRequest);

		if (areaRequest.status.poison) {

			BpCoordFileRequest fileRequest = new BpCoordFileRequest(areaRequest, null, null, areaRequest.status);
			fileRequestQueue.add(fileRequest);
			return;
		}

		fileRequestQueue.add(new BpCoordFileRequest(areaRequest, areaRequest.start, areaRequest.end, areaRequest.status));		
	}
}
