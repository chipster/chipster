package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ByteRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FileResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FsfStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class TreeNode {
	private static final long RESOLUTION = 256;
	private TreeThread tree;
	public RegionContent[] concisedValues;
	
	public BpCoord nodeBpStart;
	public BpCoord nodeBpEnd;
	public ByteRegion byteRegion;

	private TreeNode left;
	private TreeNode right;
	private TreeNode parent;
	
	private int depth;
	
	private boolean requestDistributor;

	private ByteRegion unExactByteRegion;
	private boolean isLeaf;
	private FileParser inputParser;

	public TreeNode(ByteRegion nodeByteRegion, TreeThread tree, TreeNode parent) {

		this.tree = tree;
		this.parent = parent;
		this.inputParser = tree.getInputParser();

		this.byteRegion = nodeByteRegion.clone();
		this.unExactByteRegion = nodeByteRegion.clone();

		this.isLeaf = (nodeByteRegion.getLength() <= inputParser
				.getDefaulChunkLength());
		
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

	private void fetchConcisedContent(AreaRequest areaRequest) {
		tree.createFileRequest(areaRequest, this.byteRegion, this);
	}

//	public void resetRequestDistributor() {
//		requestDistributor = false;
//		if (left != null) {
//			left.resetRequestDistributor();
//		}
//		if (right != null) {
//			right.resetRequestDistributor();
//		}
//	}

	public void processAreaRequest(AreaRequest areaRequest) {

		if (this.isLeaf) {

			if (areaRequest.status.concise) {

				createConcisedResult(areaRequest.status);
			} else {

				// Concised value isn't enough, file has to be read
				tree.createFileRequest(areaRequest, byteRegion, this);
			}
		} else {

			if ( (nodeBpStart == null || areaRequest.end.compareTo(nodeBpStart) >= 0) && 
					(nodeBpEnd == null || areaRequest.start.compareTo(nodeBpEnd) <= 0)) {

				createChildrenIfNecessary();

				if ((left.nodeBpEnd == null && right.nodeBpStart == null ) ||
						(left.nodeBpEnd != null && areaRequest.start.compareTo(left.nodeBpEnd) <= 0) ||
						(right.nodeBpStart != null && areaRequest.start.compareTo(right.nodeBpStart) > 0)) {

					if(!areaRequest.status.concise || (depth > 10 && (requestDistributor = !requestDistributor))) {
						left.processAreaRequest(areaRequest);
					}
				}

				if ((left.nodeBpEnd == null && right.nodeBpStart == null ) ||
						(left.nodeBpEnd != null && areaRequest.end.compareTo(left.nodeBpEnd) > 0) ||
						(right.nodeBpStart != null && areaRequest.end.compareTo(right.nodeBpStart) >= 0)) {

					if(!areaRequest.status.concise || (depth > 10 && !(requestDistributor = !requestDistributor))) {
						right.processAreaRequest(areaRequest);
					}
				}						
			}
		}
	}

	private void createConcisedResult(FsfStatus status) {

		for (RegionContent regCont : concisedValues) {
			tree
					.createAreaResult(new AreaResult<RegionContent>(status,
							regCont));
		}
	}

	/**
	 * @param fileResult
	 */
	public void processFileResult(FileResult fileResult) {

		if (fileResult.request.node == this) {

			if (byteRegion.exact == false) {
				byteRegion = fileResult.exactRegion;
			}

			if (concisedValues == null || nodeBpStart == null || nodeBpEnd == null) {

				FileParser parser = fileResult.chunkParser;

				nodeBpStart = parser.getBpRegion().start;
				nodeBpEnd = parser.getBpRegion().end;

				concisedValues = fileResult.chunkParser.concise(new BpCoordRegion(nodeBpStart, nodeBpEnd));

				if (fileResult.request.areaRequest != null) {
					// Continue stopped recursion now when bp location of this
					// node is known
					processAreaRequest(fileResult.request.areaRequest);
				}
			}

			if (new BpCoordRegion(nodeBpStart, nodeBpEnd).intercepts(fileResult.request.areaRequest)) {

				if (fileResult.status.concise) {
					createConcisedResult(fileResult.status);

				} else {
					getAllRows(fileResult.chunkParser,
							fileResult.request.areaRequest, fileResult.status);
				}
			}
		} else {
			updateNodeBpRegion();
		}

		if (parent != null) {
			parent.processFileResult(fileResult);
		}
	}

	private void updateNodeBpRegion() {

		if (left.nodeBpStart != null) {
			this.nodeBpStart = left.nodeBpStart;
		}
		
		if (right.nodeBpStart != null) {
			this.nodeBpStart = right.nodeBpStart;
		}				
	}

	public void getAllRows(FileParser chunkParser, AreaRequest areaRequest,
			FsfStatus status) {

		for (RegionContent rc : chunkParser
				.getAll(areaRequest.requestedContents)) {

			tree.createAreaResult(new AreaResult<RegionContent>(status, rc));
		}
	}
}
