package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.Queue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.LinkedBlockingQueue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.FastaDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;

/**
 * Processing layer of sequence data. There is no need for processing with this data type, so just passing 
 * original data through.
 * 
 * @author Petri Klemelä
 *
 */
public class FastaHandlerThread extends AreaRequestHandler {

	private FastaDataSource dataSource;
	
	private FastaFileFetcherThread fileFetcher;
	private BlockingQueue<BpCoordFileRequest> fileRequestQueue = new LinkedBlockingQueue<BpCoordFileRequest>();
	private ConcurrentLinkedQueue<ParsedFileResult> fileResultQueue = new ConcurrentLinkedQueue<ParsedFileResult>();

    public FastaHandlerThread(DataSource file, Queue<AreaRequest> areaRequestQueue,
            AreaResultListener areaResultListener) {
        
        super(areaRequestQueue, areaResultListener);
        dataSource = (FastaDataSource) file;
    }

	@Override
	public synchronized void run() {

		// Start file processing layer thread
		fileFetcher = new FastaFileFetcherThread(fileRequestQueue, fileResultQueue, this,
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
    	
    		createAreaResult(new AreaResult(fileResult.getStatus(), fileResult.getContents()));
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
