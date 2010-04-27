package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.io.IOException;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;
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

		ByteChunk chunk = new ByteChunk(inputParser.getChunkMaxByteLength());
		chunk.rowIndex = fileRequest.rowRegion.start;

		long filePosition = inputParser.getFilePosition(chunk.rowIndex);

		chunk.byteLength = dataSource.read(filePosition, chunk.byteContent);

		fileRequest.status.maybeClearQueue(fileResultQueue);
		fileRequest.status.fileRequestCount = fileRequestQueue.size();

		FileParser inputParser = (FileParser) this.inputParser.clone();
		inputParser.setChunk(chunk);

		fileResultQueue.add(new FileResult(fileRequest, inputParser, fileRequest.status));
		treeThread.notifyTree();
	}

	public long getRowCount() {
		if (this.isAlive()) {
			throw new IllegalStateException("must be called before the thread is started");
		}

		try {
			return inputParser.getRowIndex(dataSource.length());
			
		} catch (IOException e) {
			e.printStackTrace(); // FIXME fix exception handling 
		}
		return 0;
	}
}
