package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ByteRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FileResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FsfStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class TreeNode {
	private static final long RESOLUTION = 256;
	private TreeThread tree;
	public RegionContent[] concisedValues;
	public BpCoordRegion nodeBpRegion;
	public ByteRegion byteRegion;

	private TreeNode left;
	private TreeNode right;
	private TreeNode parent;

	private ByteRegion subtreeRows;
	private int maxChildCount;
	private FileParser inputParser;

	public TreeNode(ByteRegion subtreeReadIndexes, TreeThread tree, TreeNode parent) {

		this.tree = tree;
		this.parent = parent;
		this.subtreeRows = subtreeReadIndexes;
		this.inputParser = tree.getInputParser();

		this.byteRegion = inputParser.getChunkRegionMiddleOf(subtreeReadIndexes);
		this.maxChildCount = inputParser.getChildCount(subtreeReadIndexes);

		if (maxChildCount == 0) {
			this.byteRegion = subtreeReadIndexes;
		}
	}

	private void createChildren() {
		if (maxChildCount >= 1 && left == null) {
			left = new TreeNode(new ByteRegion(subtreeRows.start, byteRegion.start, false), tree, this);
		}
		if (maxChildCount == 2 && right == null) {
			right = new TreeNode(new ByteRegion(byteRegion.end, subtreeRows.end, false), tree, this);
		}
	}

	private void fetchConcisedContent(AreaRequest areaRequest) {
		tree.createFileRequest(areaRequest, this.byteRegion, this);
	}

//	 public void resetRequestDistributor() {
//		requestDistributor = false;
//		if (left != null) {
//			left.resetRequestDistributor();
//		}
//		if (right != null) {
//			right.resetRequestDistributor();
//		}
//	}

	public void processAreaRequest(AreaRequest areaRequest) {
		
		if (nodeBpRegion != null) {

			boolean recursionNeeded = areaRequest.getLength() / RESOLUTION < getSubtreeBpLengthApproximation();

			// Create file request or return value of this node
			if (nodeBpRegion.intercepts(areaRequest)) {

				if (areaRequest.status.concise) {
					createConcisedResult(areaRequest.status);
				} else {
					// Concised value isn't enough, file has to be read

					tree.createFileRequest(areaRequest, byteRegion, this);
				}
			}

			if (recursionNeeded) {
				if (maxChildCount > 0 && left == null) {
					createChildren();
				}

				if (left != null && areaRequest.start.compareTo(nodeBpRegion.start) <= 0) {
					left.processAreaRequest(areaRequest);
				}

				if (right != null && areaRequest.end.compareTo(nodeBpRegion.end) > 0) {
					right.processAreaRequest(areaRequest);
				}
			}

		} else {
			// Continue recursion when the data arrives
			fetchConcisedContent(areaRequest);
		}
	}

	private void createConcisedResult(FsfStatus status) {

		for (RegionContent regCont : concisedValues) {
			tree.createAreaResult(new AreaResult<RegionContent>(status, regCont));
		}
	}

	private long getSubtreeBpLengthApproximation() {
		if (parent != null && parent.nodeBpRegion != null && this.nodeBpRegion != null && parent.nodeBpRegion.start.chr.equals(nodeBpRegion.start.chr)) {

			return Math.abs(parent.nodeBpRegion.start.minus(this.nodeBpRegion.start) * 2);
		} else {
			return Long.MAX_VALUE;
		}
	}

	/**
	 * @param fileResult
	 */
	public void processFileResult(FileResult fileResult) {
		
		if (byteRegion.exact == false) {
			byteRegion = fileResult.exactRegion;
		}

		if (concisedValues == null || nodeBpRegion == null) {

			FileParser parser = fileResult.chunkParser;

			nodeBpRegion = parser.getBpRegion();

			concisedValues = fileResult.chunkParser.concise(nodeBpRegion);

			if (fileResult.request.areaRequest != null) {
				// Continue stopped recursion now when bp location of this node is known
				processAreaRequest(fileResult.request.areaRequest);
			}
		}

		if (nodeBpRegion.intercepts(fileResult.request.areaRequest)) {

			if (fileResult.status.concise) {
				createConcisedResult(fileResult.status);
				
			} else {
				getAllRows(fileResult.chunkParser, fileResult.request.areaRequest, fileResult.status);
			}
		}
	}

	public void getAllRows(FileParser chunkParser, AreaRequest areaRequest, FsfStatus status) {

		for (RegionContent rc : chunkParser.getAll(areaRequest.requestedContents)) {

			tree.createAreaResult(new AreaResult<RegionContent>(status, rc));
		}
	}
}
