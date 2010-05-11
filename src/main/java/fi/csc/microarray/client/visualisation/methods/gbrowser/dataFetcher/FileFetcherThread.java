package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.io.IOException;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ByteRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FileRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FileResult;

public class FileFetcherThread extends Thread {

	private BlockingQueue<FileRequest> fileRequestQueue;
	private ConcurrentLinkedQueue<FileResult> fileResultQueue;

	private TreeThread treeThread;

	private DataSource dataSource;

	private FileParser inputParser;

	public FileFetcherThread(BlockingQueue<FileRequest> fileRequestQueue, ConcurrentLinkedQueue<FileResult> fileResultQueue, TreeThread treeThread, FileParser inputParser) {

		this.fileRequestQueue = fileRequestQueue;
		this.fileResultQueue = fileResultQueue;
		this.treeThread = treeThread;
		this.inputParser = inputParser;

		this.setDaemon(true);

		this.dataSource = treeThread.getFile();
	}

	public void run() {

		while (true) {
			try {
				processFileRequest(fileRequestQueue.take());
				// if(fileRequestQueue.peek() != null){
				// processFileRequest(fileRequestQueue.poll());
				// }
			} catch (IOException e) {
				e.printStackTrace(); // FIXME fix exception handling
			} catch (InterruptedException e) {
				e.printStackTrace(); // FIXME fix exception handling
			}
		}
	}
	
	private void processFileRequest(FileRequest fileRequest) throws IOException {
		
		String chunk;
		ByteRegion exactRegion = null;
		
		if (fileRequest.byteRegion.exact) {

			byte[] byteChunk = new byte[(int)fileRequest.byteRegion.getLength()];
				
			dataSource.read(fileRequest.byteRegion.start, byteChunk);			
			
			chunk = new String(byteChunk);

		} else {
			
			//some extra to get the last line fully
			byte[] byteChunk = new byte[(int)fileRequest.byteRegion.getLength() + 1000];
			
			int length = dataSource.read(fileRequest.byteRegion.start, byteChunk);
			
			String file = new String(byteChunk).substring(0, length);		
			
			exactRegion = new ByteRegion();
			int i = 0;
			
			if(fileRequest.byteRegion.start != 0) {
				for (; ; i++) {
					
					if ( i >= file.length()) {
						//not a single new line found, source file is broken
						return;
					}
					
					if (file.charAt(i) == '\n') {
						i++;
						exactRegion.start = fileRequest.byteRegion.start + i;
						break;
					}
				}
			}
			
			StringBuffer lines = new StringBuffer();
			
			for (; ; i++) {
				
				
				lines.append(file.charAt(i));		
				
				if (file.charAt(i) == '\n' && i >= fileRequest.byteRegion.getLength()) {
					break;
				}
				
				if ( i >= file.length() - 1) {	
					
					//TODO buffer ended before the new line chracter, discard the last line 
					return;
				}
			}
			
			exactRegion.end = fileRequest.byteRegion.start + i;			
			exactRegion.exact = true;

			chunk = lines.toString();				
		}

		fileRequest.status.maybeClearQueue(fileResultQueue);
		fileRequest.status.fileRequestCount = fileRequestQueue.size();

		FileResult result = new FileResult(chunk, fileRequest, inputParser, exactRegion, fileRequest.status);

		fileResultQueue.add(result);
		treeThread.notifyTree();

	}

	public long getFileLength() {
		if (this.isAlive()) {
			throw new IllegalStateException("must be called before the thread is started");
		}

		try {
			return dataSource.length();
		} catch (IOException e) {
			e.printStackTrace(); // FIXME fix exception handling 
		}
		return 0;
	}
}
