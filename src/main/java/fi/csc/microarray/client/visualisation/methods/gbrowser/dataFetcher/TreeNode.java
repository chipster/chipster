package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FileResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FsfStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RowRegion;

public class TreeNode {
	private static final long RESOLUTION = 256;
	private TreeThread tree;
	public RegionContent[] concisedValues;
	public BpCoordRegion nodeBpRegion;
	public RowRegion rowRegion;

	private TreeNode left;
	private TreeNode right;
	private TreeNode parent;

	private RowRegion subtreeRows;
	private int maxChildCount;
	private FileParser inputParser;

	public TreeNode(RowRegion subtreeReadIndexes, TreeThread tree, TreeNode parent) {

		this.tree = tree;
		this.parent = parent;
		this.subtreeRows = subtreeReadIndexes;
		this.inputParser = tree.getInputParser();

		this.rowRegion = inputParser.getChunkRegionMiddleOf(subtreeReadIndexes);
		this.maxChildCount = inputParser.getChildCount(subtreeReadIndexes);

		if (maxChildCount == 0) {
			this.rowRegion = subtreeReadIndexes;
		}
	}

	private void createChildren() {
		if (maxChildCount >= 1 && left == null) {
			left = new TreeNode(new RowRegion(subtreeRows.start, rowRegion.start - 1), tree, this);
		}
		if (maxChildCount == 2 && right == null) {
			right = new TreeNode(new RowRegion(rowRegion.end + 1, subtreeRows.end), tree, this);
		}
	}

	private void fetchConcisedContent(AreaRequest areaRequest) {
		tree.createFileRequest(areaRequest, this.rowRegion, this);
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

					tree.createFileRequest(areaRequest, rowRegion, this);
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
			// Continue recursion when to data arrives
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
	 *            null if should be calculated from children
	 */
	public void processFileResult(FileResult fileResult) {

		if (concisedValues == null || nodeBpRegion == null) {

			FileParser parser = fileResult.chunkParser;
			Long row = parser.getRowIndex();

			Chromosome chr = (Chromosome) parser.get(row, ColumnType.CHROMOSOME);

			BpCoord start = new BpCoord((Long) parser.get(row, ColumnType.BP_START), chr);
			BpCoord end = new BpCoord((Long) parser.get(row + parser.getChunkRowCount() - 1, ColumnType.BP_START), chr);

			nodeBpRegion = new BpCoordRegion(start, end);

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

		for (int i = 0; i < chunkParser.getChunkRowCount(); i++) {

			BpCoordRegion bpRegion = chunkParser.getBpRegion(i + chunkParser.getRowIndex());

			Map<ColumnType, Object> values = chunkParser.getValues(i + chunkParser.getRowIndex(), areaRequest.requestedContents);

			values.put(ColumnType.FILE_INDEX, i + chunkParser.getRowIndex());
			tree.createAreaResult(new AreaResult<RegionContent>(status, new RegionContent(bpRegion, values)));
		}
	}
}
