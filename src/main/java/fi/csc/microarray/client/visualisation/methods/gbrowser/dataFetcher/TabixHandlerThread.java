package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.Queue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.LinkedBlockingQueue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ParsedFileResult;


public abstract class TabixHandlerThread extends AreaRequestHandler {
    
	protected TabixDataSource dataSource;
	protected BlockingQueue<BpCoordFileRequest> fileRequestQueue = new LinkedBlockingQueue<BpCoordFileRequest>();
	protected ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue = new ConcurrentLinkedQueue<ParsedFileResult>();

    public TabixHandlerThread(DataSource file, Queue<AreaRequest> areaRequestQueue,
            AreaResultListener areaResultListener) {
        
        super(areaRequestQueue, areaResultListener);
        dataSource = (TabixDataSource) file;
    }

//	@Override
//	public synchronized void run() {
//
//		// Start file processing layer thread
//		fileFetcher = new BedTabixFileFetcherThread(fileRequestQueue, fileResultQueue, this, tabixData);
//		fileFetcher.start();
//		
//		// Start this thread
//		super.run();
//	}

	protected boolean checkOtherQueues() {
		ParsedFileResult fileResult = null;
		if ((fileResult = fileResultQueue.poll()) != null) {
			processFileResult(fileResult);
		}
		return fileResult != null;
	}

	protected void processFileResult(ParsedFileResult fileResult) {

		createAreaResult(new AreaResult(fileResult.getStatus(), fileResult.getContents()));
	}

    @Override
    protected void processAreaRequest(AreaRequest areaRequest) {
    	
		super.processAreaRequest(areaRequest);
		
		if (areaRequest.getStatus().poison) {
			
			BpCoordFileRequest fileRequest = new BpCoordFileRequest(areaRequest, null, null, areaRequest.getStatus());
			fileRequestQueue.add(fileRequest);
			return;
		}
    	
		fileRequestQueue.add(new BpCoordFileRequest(areaRequest, areaRequest.start, areaRequest.end, areaRequest.getStatus()));
    }
}
