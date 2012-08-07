package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.Queue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.LinkedBlockingQueue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.IndexedFastaDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;



public class IndexedFastaHandlerThread extends AreaRequestHandler {
    
	private IndexedFastaDataSource data;
	private IndexedFastaFileFetcherThread fileFetcher;
	private BlockingQueue<BpCoordFileRequest> fileRequestQueue = new LinkedBlockingQueue<BpCoordFileRequest>();
	private ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue = new ConcurrentLinkedQueue<ParsedFileResult>();

    public IndexedFastaHandlerThread(DataSource file, Queue<AreaRequest> areaRequestQueue,
            AreaResultListener areaResultListener) {
        
        super(areaRequestQueue, areaResultListener);
        data = (IndexedFastaDataSource) file;
    }

	@Override
	public synchronized void run() {

		// Start file processing layer thread
		fileFetcher = new IndexedFastaFileFetcherThread(fileRequestQueue, fileResultQueue, this, data);
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

      		createAreaResult(new AreaResult(fileResult.getStatus(), fileResult.getContents()));
	}

	/**
     * Handles normal and concised area requests by using SAMFile.
     */
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
