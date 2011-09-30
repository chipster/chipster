package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ByteRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ChunkFileResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FsfStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class TreeNode {

	private static final boolean DEPTH_LIMIT_ACTIVE = false;
	
	private ChunkTreeHandlerThread tree;
	public RegionContent[] concisedValues;

	public BpCoord nodeBpStart;
	public ByteRegion byteRegion;

	private TreeNode left;
	private TreeNode right;
	private TreeNode parent;

	private int depth;

	private ByteRegion unExactByteRegion;
	private boolean isLeaf;

	/**
	 * DOCME
	 * 
	 * Chunk reader is the original genome browser data fetching implementation with a tree type 
	 * data structure and it was made primarily for the tab separated text files. This class is the
	 * implementation of a node of tree structure. All the recursive functionality is done in this 
	 * class. Things concerning the whole tree are done in class called ChunkTreeHandlerThread.
	 * 
	 * Only the leaf nodes of the tree point to the part of the file and the upper nodes of the tree
	 * are used only for the searching and sampling.
	 * 
	 */
	public TreeNode(ByteRegion nodeByteRegion, ChunkTreeHandlerThread tree, TreeNode parent) {

		this.tree = tree;
		this.parent = parent;

		this.byteRegion = nodeByteRegion.clone();
		this.unExactByteRegion = nodeByteRegion.clone();

		this.isLeaf = (nodeByteRegion.getLength() <= tree.getInputParser().getDefaulChunkLength());

		//The depth isn't really needed, but makes debugging easier
		if (parent == null) {
			depth = 0;
		} else {
			depth = parent.depth + 1;
		}
	}

	/**
	 * If this is leaf and at least one children is null, create children and split the area for them.
	 */
	private void createChildrenIfNecessary() {
		if (!isLeaf && (left == null || right == null)) {
			left = new TreeNode(new ByteRegion(unExactByteRegion.start, (long) unExactByteRegion.getMid() - 1, false), tree, this);
			right = new TreeNode(new ByteRegion((long) unExactByteRegion.getMid(), unExactByteRegion.end, false), tree, this);
		}
	}

	/**
	 * Recursive method to find out the genomic location of the most left leaf under this node.
	 * 
	 * @param areaRequest
	 * @param source
	 */
	private void updateNodeBpStart(AreaRequest areaRequest, TreeNode source) {
		if (this.isLeaf) {
			areaRequest.status.bpSearchSource = source;
			tree.createFileRequest(areaRequest, this.byteRegion, this);
			
		} else {
			createChildrenIfNecessary();
			left.updateNodeBpStart(areaRequest, source);
		}
	}

	/**
	 * Methods ProcessAreaRequest and ProcessFileResult implement the most essential part of the
	 * technical performance optimizations of the genome browser, called dynamic search tree and 
	 * gradual sampling, see Chipster wiki for more information about these. Understanding these 
	 * main concepts is essential to interpret the complicated interaction between
	 * these methods.
	 * 
	 * @param areaRequest
	 */
	public void processAreaRequest(AreaRequest areaRequest) {

		// if on leaf, do not recurse down but read file (if needed) and return result
		if (this.isLeaf) {

			// concised data requested
			if (areaRequest.status.concise) {

				if (concisedValues == null) {
					// create file request to get the summary
					tree.createFileRequest(areaRequest, this.byteRegion, this);
					
				} else {
					// we have the summary already, return it
					createConcisedResult(areaRequest, areaRequest.status);
				}
				
			} else {
				// non-concised result wanted
				tree.createFileRequest(areaRequest, byteRegion, this);
			}
		} else {

			// create childrens of the node if they aren't there already
			createChildrenIfNecessary();


			/* The genomic location of the most left leaf under the right child is the dividing 
			 * point between the left and right children and needed for the searching. If we don't
			 * have this data yet, we have to stop execution here and continue later when the 
			 * FileResult arrives.
			 */ 
			if (right.nodeBpStart == null) {
				right.updateNodeBpStart(areaRequest, this);

			} else {
				// recurse down

				try {
					
					/* limit splitting of sampling to certain count, effectively this limits the
					 * amount of sampling results. Sampling is using the concised data. If the 
					 * concised data isn't enough, we can't limit the searching.
					 */	
					boolean canSplit;
					if (areaRequest.status.concise) {
						canSplit = !DEPTH_LIMIT_ACTIVE || areaRequest.depthToGo > 0;
					} else {
						canSplit = true;
					}
					
					boolean recurseLeft = areaRequest.start.compareTo(right.nodeBpStart) < 0; 
					boolean recurseRight = areaRequest.end.compareTo(right.nodeBpStart) > 0;
					boolean recurseBoth = recurseLeft && recurseRight;
					
					// solve conflicts with dice
					if (!canSplit && recurseBoth) {
						// can't recurse to both directions because splitting forbidden
						if (Math.random() < 0.5d) {
							recurseLeft = false;
						} else {
							recurseRight = false;
						}
					}
					
					// recurse to left
					if (recurseLeft) {
						AreaRequest clone = areaRequest.clone();
						if (recurseBoth) {
							clone.depthToGo--;
						}
						left.processAreaRequest(clone);
					}

					// recurse to right
					if (recurseRight) {
						AreaRequest clone = areaRequest.clone();
						if (recurseBoth) {
							clone.depthToGo--;
						}
						right.processAreaRequest(clone);
					}
					
				} catch (CloneNotSupportedException e) {
					throw new RuntimeException(e);
				}
			}
		}
	}

	/**
	 * Methods ProcessAreaRequest and ProcessFileResult implement the most essential part of the
	 * technical performance optimizations of the genome browser, called dynamic search tree and 
	 * gradual sampling, see Chipster wiki for more information about these. Understanding these 
	 * main concepts is essential to interpret the complicated interaction between
	 * these methods.
	 * 
	 * @See FsfStatus for description of chunk
	 * @param fileResult
	 */
	public void processFileResult(ChunkFileResult fileResult) {

		if (isLeaf) {

			/* Store the location of full lines now when it's known, see ChunkFileFetcherThread
			 * for more information about this optimization.
			 */			
			if (byteRegion.exact && fileResult.exactRegion != null) {
				byteRegion = fileResult.exactRegion;
			}

			/* If the data of this leaf is missing, save it now to make further searching and 
			 * concised requests quicker.
			 */
			if (concisedValues == null || nodeBpStart == null) {

				FileParser parser = fileResult.chunkParser;

				nodeBpStart = parser.getBpRegion(fileResult.chunk).start;

				if (parent != null) {
					parent.nodeBpStartUpdated(this);
				}

				concisedValues = fileResult.chunkParser.concise(fileResult.chunk);
			}

			/* Create the result object if the result intercepts with the requested area. 
			 * Normally this is always true, but it's not a big job to check anyway and make it
			 * little bit more robust.
			 */
			if (fileResult.request.areaRequest.intersects(fileResult.chunkParser.getBpRegion(fileResult.chunk))) {

				if (fileResult.status.concise) {
					createConcisedResult(fileResult.request.areaRequest, fileResult.status);

				} else {
					createAreaResultOfAllRows(fileResult.chunk, fileResult.chunkParser, fileResult.request.areaRequest, fileResult.status);
				}
			}

		} else {

			/* If the genomic location of node wasn't known, it was requested from the 
			 * updateNodeBpStart method. In this case, continue searching.
			 */
			if (fileResult.status.bpSearchSource == this) {

				fileResult.status.bpSearchSource = null;

				// Continue finding of specific place in file now when the location of this branch is known
				processAreaRequest(fileResult.request.areaRequest);
			}
		}

		if (parent != null) {
			//TODO might be useless after the bpSeachSource has been found
			parent.processFileResult(fileResult);
		}
	}

	/**
	 * When the genomic location of node is found, for example after a call of updateNodeBpStart 
	 * method, it will be updated upwards with this method if it's the starting position for
	 * other nodes also.
	 * 
	 * @param source
	 */
	private void nodeBpStartUpdated(TreeNode source) {
		// update start coordinate and recurse up the tree
		if (source == left) {
			this.nodeBpStart = left.nodeBpStart;
			if (parent != null) {
				parent.nodeBpStartUpdated(this);
			}
		}
	}

	/**
	 * Create concised results from data that we had in a memory or we just read from the file.
	 * 
	 * @param areaRequest
	 * @param status
	 */
	private void createConcisedResult(AreaRequest areaRequest, FsfStatus status) {

		LinkedList<RegionContent> contents = new LinkedList<RegionContent>();
		
		for (RegionContent regCont : concisedValues) {
			if (areaRequest.intersects(regCont.region)) {
				contents.add(regCont);
			}
		}
		
		if (!contents.isEmpty()) {
			tree.createAreaResult(new AreaResult(status, contents));
		}

	}

	/**
	 * Create the result containing all the rows from the file intercepting the request area.
	 * 
	 * @param chunk
	 * @param chunkParser
	 * @param areaRequest
	 * @param status
	 */
	public void createAreaResultOfAllRows(Chunk chunk, FileParser chunkParser, AreaRequest areaRequest, FsfStatus status) {

		LinkedList<RegionContent> contents = new LinkedList<RegionContent>();

		for (RegionContent rc : chunkParser.getAll(chunk, areaRequest.requestedContents)) {
			if (areaRequest.intersects(rc.region)) {
				contents.add(rc);
			}
		}

		if (!contents.isEmpty()) {
			tree.createAreaResult(new AreaResult(status, contents));
		}
	}
}
