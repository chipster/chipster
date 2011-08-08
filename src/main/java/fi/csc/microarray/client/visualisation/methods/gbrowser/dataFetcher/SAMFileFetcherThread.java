package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.io.IOException;
import java.util.List;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.SAMDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;


public class SAMFileFetcherThread extends Thread {

	private BlockingQueue<SAMFileRequest> fileRequestQueue;
	private ConcurrentLinkedQueue<SAMFileResult> fileResultQueue;

	private SAMDataSource dataSource;

	public SAMFileFetcherThread(BlockingQueue<SAMFileRequest> fileRequestQueue,
	        ConcurrentLinkedQueue<SAMFileResult> fileResultQueue,
	        SAMDataSource dataSource) {

		this.fileRequestQueue = fileRequestQueue;
		this.fileResultQueue = fileResultQueue;
		this.dataSource = dataSource;
		this.setDaemon(true);
	}

	public void run() {

		while (true) {
			try {
				processFileRequest(fileRequestQueue.take());
				
			} catch (IOException e) {
				e.printStackTrace(); // FIXME fix exception handling
			} catch (InterruptedException e) {
				e.printStackTrace(); // FIXME fix exception handling
			}
		}
	}

	private void processFileRequest(SAMFileRequest fileRequest) throws IOException {
		List<RegionContent> content;
		if (fileRequest.areaRequest.status.concise) { 
			content = dataSource.getSAM().getConciseReads(fileRequest.areaRequest);
			
		} else {
			content = dataSource.getSAM().getReads(fileRequest.areaRequest);
		}
		
		SAMFileResult result = new SAMFileResult(content, fileRequest, fileRequest.areaRequest, fileRequest.getStatus());

		fileResultQueue.add(result);
		
		synchronized (this) {
			notifyAll();
		}
		
	}
}
