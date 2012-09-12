package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.io.IOException;
import java.util.Queue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.LinkedBlockingQueue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ChunkFileRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ChunkFileResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ByteRegion;

/**
 * Processing layer implementation for chunk files (tab separated text files). Because those files do not 
 * have indexing, we have to build a tree index for them in the fly. The tree does two things: it builds
 * an index that can map genomic coordinates to file coordinates, and summarises values so that high level
 * views can be drawn by using the tree alone.
 * 
 * @author Petri Klemelä
 * @see TreeNode
 */
public class ChunkTreeHandlerThread extends AreaRequestHandler {

	private TreeNode rootNode;

	private BlockingQueue<ChunkFileRequest> fileRequestQueue = new LinkedBlockingQueue<ChunkFileRequest>();
	private ConcurrentLinkedQueue<ChunkFileResult> fileResultQueue = new ConcurrentLinkedQueue<ChunkFileResult>();

	private ChunkFileFetcherThread fileFetcher;

	private static final boolean DEBUG = false;

	private FileParser inputParser;

	private ChunkDataSource file;

	public ChunkTreeHandlerThread(DataSource file, Queue<AreaRequest> areaRequestQueue,
	        AreaResultListener areaResultListener) {
	    
		super(areaRequestQueue, areaResultListener);
		this.file = (ChunkDataSource) file;
		this.inputParser = this.file.getFileParser();
	}

	@Override
	public synchronized void run() {

		try {
			fileFetcher = new ChunkFileFetcherThread(fileRequestQueue, fileResultQueue, this, inputParser);
			createTree(fileFetcher.getFileLength());
			fileFetcher.start();

			super.run();

		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}

	protected boolean checkOtherQueues() {
		ChunkFileResult fileResult = null;
		if ((fileResult = fileResultQueue.poll()) != null) {
			fileResult.status.fileResultCount = fileResultQueue.size();
			processFileResult(fileResult);
		}
		return fileResult != null;
	}

	private void createTree(long fileLength) throws IOException {
		rootNode = new TreeNode(new ByteRegion(file.getHeaderLength(), fileLength, false), this, null);
	}

	private void processFileResult(ChunkFileResult fileResult) {
		fileResult.request.node.processFileResult(fileResult);
	}

	protected void processAreaRequest(AreaRequest areaRequest) {
		
		super.processAreaRequest(areaRequest);
		
		if (areaRequest.status.poison) {
			
			ChunkFileRequest fileRequest = new ChunkFileRequest(areaRequest, null, null, areaRequest.status);
			fileRequestQueue.add(fileRequest);
			return;
		}
		
		if (DEBUG) {
			System.out.println("Tree: Got area request " + areaRequest);
		}
		rootNode.processAreaRequest(areaRequest);
	}

	/**
	 * @param areaRequest
	 * @param byteRegion will be cloned to allow modifications in other thread without causing 
	 * problems here
	 * @param node
	 */
	public void createFileRequest(AreaRequest areaRequest, ByteRegion byteRegion, TreeNode node) {
		areaRequest.status.maybeClearQueue(fileRequestQueue);

		fileRequestQueue.add(new ChunkFileRequest(areaRequest, byteRegion.clone(), node, areaRequest.status));
	}

	public FileParser getInputParser() {
		return inputParser;
	}

	public ChunkDataSource getFile() {
		return file;
	}
	
	public String toString() {
		return this.getClass().getName() + " - " + file;
	}
}
