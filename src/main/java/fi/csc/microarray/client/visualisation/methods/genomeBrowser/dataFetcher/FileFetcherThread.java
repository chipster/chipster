package fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.ReadInstructions;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.FileRequest;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.FileResult;


public class FileFetcherThread extends Thread{
				
	private BlockingQueue<FileRequest> fileRequestQueue;
	private ConcurrentLinkedQueue<FileResult> fileResultQueue;
	
	private TreeThread treeThread;
	
	private RandomAccessFile raf;
	private ReadInstructions<?> instructions;

	public FileFetcherThread(BlockingQueue<FileRequest> fileRequestQueue,
			ConcurrentLinkedQueue<FileResult> fileResultQueue, TreeThread treeThread, ReadInstructions<?> instructions) {
		this.fileRequestQueue = fileRequestQueue;
		this.fileResultQueue = fileResultQueue;
		this.treeThread = treeThread;
		this.instructions = instructions;
		
		this.setDaemon(true);
		
		try {
			raf = new RandomAccessFile(treeThread.getFile(), "r");
			//raf = new RandomAccessFile(new File("bowtie.fsf"), "r");
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		}
	}

	public void run() {
		while(true){
			try {
				
				processFileRequest(fileRequestQueue.take());
				
//				if(fileRequestQueue.peek() != null){
//					processFileRequest(fileRequestQueue.poll());					
//				}
			} catch (IOException e) {
				e.printStackTrace();
			} catch (InterruptedException e) {
					e.printStackTrace();
			}						
		}
	}

	private void processFileRequest(FileRequest fileRequest) throws IOException {
		
//		System.out.println("File: Got file request " + fileRequest.region.start);
		
		ByteChunk chunk = new ByteChunk(instructions.chunker.getChunkByteLength());			
		chunk.readIndex = fileRequest.rowRegion.start;

		raf.seek(instructions.chunker.getFilePosition(chunk.readIndex));						
		chunk.length = raf.read(chunk.content);
		
		fileRequest.status.maybeClearQueue(fileResultQueue);
		fileRequest.status.fileRequestCount = fileRequestQueue.size();
		
		fileResultQueue.add(new FileResult(fileRequest, chunk, fileRequest.status));
		treeThread.notifyTree();
		
	}

	public long getRowCount() {
		if(this.isAlive()){
			throw new IllegalStateException("Must be called before the thread is started");
		}
		
		try {
			return instructions.chunker.getReadIndex(raf.length());
		} catch (IOException e) {
			e.printStackTrace();
		}
		return 0;
	}

}
