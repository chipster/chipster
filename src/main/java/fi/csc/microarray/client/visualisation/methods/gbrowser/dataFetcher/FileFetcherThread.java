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
			
			dataSource.seek(fileRequest.byteRegion.start);

			byte[] byteChunk = new byte[(int)fileRequest.byteRegion.getLength()];
				
			dataSource.read(byteChunk);			
			
			chunk = new String(byteChunk);
			
			dataSource.clean();

		} else {
			
			dataSource.seek(inputParser.getFilePosition(fileRequest.byteRegion.start));
			
			exactRegion = new ByteRegion();
			
			//Find next new line
			if(fileRequest.byteRegion.start != 0) {
				dataSource.readLine();
			}

			exactRegion.start = dataSource.getPosition();
			
			StringBuilder lines = new StringBuilder();

			while (dataSource.getPosition() <= fileRequest.byteRegion.end) {

				lines.append(dataSource.readLine()); // FIXME out of mem
				lines.append("\n");
			}

			exactRegion.end = dataSource.getPosition() - 1;
			exactRegion.exact = true;

			chunk = lines.toString();				
			
			dataSource.clean();
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
