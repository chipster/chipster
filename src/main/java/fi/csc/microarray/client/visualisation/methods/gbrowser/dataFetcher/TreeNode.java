package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ByteRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FileResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FsfStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class TreeNode {

	private TreeThread tree;
	public RegionContent[] concisedValues;
	
	public BpCoord nodeBpStart;
	public ByteRegion byteRegion;

	private TreeNode left;
	private TreeNode right;
	private TreeNode parent;
	
	private int depth;
	
	private boolean requestDistributor;

	private ByteRegion unExactByteRegion;
	private boolean isLeaf;

	public TreeNode(ByteRegion nodeByteRegion, TreeThread tree, TreeNode parent) {

		this.tree = tree;
		this.parent = parent;

		this.byteRegion = nodeByteRegion.clone();
		this.unExactByteRegion = nodeByteRegion.clone();

		this.isLeaf = (nodeByteRegion.getLength() <= tree.getInputParser().getDefaulChunkLength());
		
		if (parent == null) {
			depth = 0;
		} else {
			depth = parent.depth + 1;
		}
	}

	private void createChildrenIfNecessary() {
		if (!isLeaf && (left == null || right == null)) {

			left = new TreeNode(new ByteRegion(unExactByteRegion.start,
					(long) unExactByteRegion.getMid() - 1, false), tree, this);
			right = new TreeNode(new ByteRegion((long) unExactByteRegion
					.getMid(), unExactByteRegion.end, false), tree, this);
		}
	}
	
	private void updateNodeBpStart(AreaRequest areaRequest, TreeNode source) {
		if(this.isLeaf) {
			
			areaRequest.status.bpSearchSource = source;
			tree.createFileRequest(areaRequest, this.byteRegion, this);
		} else {
			
			createChildrenIfNecessary();
			left.updateNodeBpStart(areaRequest, source);			
		}
	}

	public void processAreaRequest(AreaRequest areaRequest) {

		if (this.isLeaf) {

			if (areaRequest.status.concise) {

				if(concisedValues == null) {
					tree.createFileRequest(areaRequest, this.byteRegion, this);			
				} else {
					
					createConcisedResult(areaRequest, areaRequest.status);
				}
			} else {

				// Concised value isn't enough, file has to be read
				tree.createFileRequest(areaRequest, byteRegion, this);
			}
		} else {

			createChildrenIfNecessary();
			
			if (right.nodeBpStart == null) {
				
				right.updateNodeBpStart(areaRequest, this);
				
			} else {
									
				if( areaRequest.start.compareTo(right.nodeBpStart) < 0 &&
						(!areaRequest.status.concise || (depth < 10 || (requestDistributor = !requestDistributor)))) {

					left.processAreaRequest(areaRequest);
				}

				if( areaRequest.end.compareTo(right.nodeBpStart) > 0 &&
						(!areaRequest.status.concise || (depth < 10 || !(requestDistributor = !requestDistributor)))) {

					right.processAreaRequest(areaRequest);
				}					
			}
		}
	}



	/**
	 * @param fileResult
	 */
	public void processFileResult(FileResult fileResult) {
		
		if (isLeaf) {

			if (byteRegion.exact && fileResult.exactRegion != null) {
				byteRegion = fileResult.exactRegion;
			}

			if (concisedValues == null || nodeBpStart == null ) {

				FileParser parser = fileResult.chunkParser;

				nodeBpStart = parser.getBpRegion(fileResult.chunk).start;
				
				if (parent != null) {
					parent.updateNodeBpStart(this);
				}

				concisedValues = fileResult.chunkParser.concise(fileResult.chunk);
			}

			if (fileResult.request.areaRequest.intercepts(
					fileResult.chunkParser.getBpRegion(fileResult.chunk))) {

				if (fileResult.status.concise) {
					createConcisedResult(fileResult.request.areaRequest, fileResult.status);

				} else {
					createAreaResultOfAllRows(fileResult.chunk, fileResult.chunkParser,
							fileResult.request.areaRequest, fileResult.status);
				}
			} 
			
		} else {
			
			if (fileResult.status.bpSearchSource == this) {
				
				//fileResult.status.bpSearchSource = null;
				
				//Continue finding of specific place in file now when the location of this branch is known
				processAreaRequest(fileResult.request.areaRequest);
			}
		}

		if (parent != null) {
			parent.processFileResult(fileResult);
		}
	}
	
	private void updateNodeBpStart(TreeNode source) {
		if (source == left) {
			this.nodeBpStart = left.nodeBpStart;
			if (parent != null) {
				parent.updateNodeBpStart(this);
			}
		}
	}
	
	private void createConcisedResult(AreaRequest areaRequest, FsfStatus status) {

		for (RegionContent regCont : concisedValues) {
			
			if (areaRequest.intercepts(regCont.region)) {
				tree.createAreaResult(new AreaResult<RegionContent>(status, regCont));
			}
		}
	}

	public void createAreaResultOfAllRows(String chunk, FileParser chunkParser, AreaRequest areaRequest,
			FsfStatus status) {		

		for (RegionContent rc : chunkParser.getAll(chunk, areaRequest.requestedContents)) {
			
			if (areaRequest.intercepts(rc.region)) {
				tree.createAreaResult(new AreaResult<RegionContent>(status, rc));
			} 
		}				
	}
}
