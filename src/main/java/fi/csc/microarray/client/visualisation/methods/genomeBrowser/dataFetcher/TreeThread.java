package fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher;
import java.io.File;
import java.util.Queue;
import java.util.concurrent.BlockingQueue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.LinkedBlockingQueue;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.ReadInstructions;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.FileRequest;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.FileResult;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.FsfStatus;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.Region;


public class TreeThread<T> extends AreaRequestHandler {

	private TreeNode<T> rootNode;
	
	private BlockingQueue<FileRequest> fileRequestQueue = new LinkedBlockingQueue<FileRequest>();
	private ConcurrentLinkedQueue<FileResult> fileResultQueue = new ConcurrentLinkedQueue<FileResult>();
		
	private FileFetcherThread fileFetcher;
	
	private boolean debug = false;

	private ReadInstructions<T> instructions;
	

	private File file;
	
	public TreeThread(
			File file, 
			Queue<AreaRequest> areaRequestQueue,
			AreaResultListener areaResultListener, 
			ReadInstructions<T> instructions) {
		
		super(areaRequestQueue, areaResultListener);
		this.instructions = instructions;
		this.file = file;
	}

	
	public synchronized void run(){
		fileFetcher = new FileFetcherThread(fileRequestQueue, fileResultQueue, this, instructions);
		createTree(fileFetcher.getRowCount());
		fileFetcher.start();
		
		super.run();
	}
	
	protected boolean checkOtherQueues() {
		FileResult fileResult;
		if((fileResult = fileResultQueue.poll()) != null){
			fileResult.status.fileResultCount = fileResultQueue.size();
			
			processFileResult(fileResult);
		}
		return fileResult != null;
	}
	
	private void createTree(long rowCount) {
	
		rootNode = new TreeNode<T>(new Region(0, rowCount), this, null);
	}

	private void processFileResult(FileResult fileResult) {
		fileResult.request.node.processFileResult(fileResult);
	}

	protected void processAreaRequest(AreaRequest areaRequest) {
		if(debug)System.out.println("Tree: Got area request " + areaRequest);
		
		rootNode.processAreaRequest(areaRequest);
	}

	public void createFileRequest(AreaRequest areaRequest, Region rowRegion, TreeNode node, FsfStatus status) {
		if(debug)System.out.println("Tree: Creating file request");
		status.maybeClearQueue(fileRequestQueue);
				
		fileRequestQueue.add(new FileRequest(areaRequest, rowRegion, node, status));
	}
	
	public ReadInstructions<T> getInstructions() {
		return instructions;
	}
	
	
	public File getFile() {
		return file;
	}
}
