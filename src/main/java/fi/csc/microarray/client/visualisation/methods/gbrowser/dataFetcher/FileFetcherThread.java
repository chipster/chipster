package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ByteRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FileRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FileResult;

public class FileFetcherThread extends Thread {

	private BlockingQueue<FileRequest> fileRequestQueue;
	private ConcurrentLinkedQueue<FileResult> fileResultQueue;

	private TreeThread treeThread;

	private RandomAccessFile raf;

	private FileParser inputParser;

	public FileFetcherThread(BlockingQueue<FileRequest> fileRequestQueue, ConcurrentLinkedQueue<FileResult> fileResultQueue, TreeThread treeThread, FileParser inputParser) {

		this.fileRequestQueue = fileRequestQueue;
		this.fileResultQueue = fileResultQueue;
		this.treeThread = treeThread;
		this.inputParser = inputParser;

		this.setDaemon(true);

		try {
			raf = new RandomAccessFile(treeThread.getFile(), "r");

		} catch (FileNotFoundException e) {
			e.printStackTrace(); // FIXME fix exception handling
		}
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

		raf.seek(inputParser.getFilePosition(fileRequest.byteRegion.start));

		ByteRegion exactRegion = new ByteRegion();

		if (fileRequest.byteRegion.exact) {

			exactRegion.start = raf.getFilePointer();

			byte[] byteChunk = new byte[(int)fileRequest.byteRegion.getLength()];
							
			raf.read(byteChunk);
			
			chunk = new String(byteChunk);

			exactRegion = fileRequest.byteRegion;

		} else {
			
			//Find next new line
			if(fileRequest.byteRegion.start != 0) {
				raf.readLine();
			}

			exactRegion.start = raf.getFilePointer();
			
			StringBuilder lines = new StringBuilder();

			while (raf.getFilePointer() <= fileRequest.byteRegion.end) {

				lines.append(raf.readLine());
				lines.append("\n");
			}

			exactRegion.end = raf.getFilePointer() - 1;
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
			return raf.length();

		} catch (IOException e) {
			e.printStackTrace(); // FIXME fix exception handling 
		}
		return 0;
	}
}
