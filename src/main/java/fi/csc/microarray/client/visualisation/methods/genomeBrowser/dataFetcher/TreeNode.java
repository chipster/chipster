package fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.ChunkParser;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.Content;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.ReadInstructions;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.FileResult;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.FsfStatus;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.RegionContent;


public class TreeNode<T>{
	private static final long RESOLUTION = 256;
	private TreeThread<T> tree;
	public T concisedValue;
	public Region nodeBpRegion;
	public Region nodeRows;

	private TreeNode<T> left;
	private TreeNode<T> right;
	private TreeNode<T> parent;

	private boolean requestDistributor = false;
	
	private Region subtreeRows;
	private int maxChildCount;
	private ReadInstructions<T> instructions;

	public TreeNode(Region subtreeReadIndexes, TreeThread<T> tree, TreeNode<T> parent) {

		this.tree = tree;
		this.parent = parent;		
		this.subtreeRows = subtreeReadIndexes;
		this.instructions = tree.getInstructions();

		this.nodeRows = new Region();

		this.nodeRows = instructions.chunker.getChunkMiddleOf(subtreeReadIndexes);
		this.maxChildCount = instructions.chunker.getChildCount(subtreeReadIndexes);

		if(maxChildCount == 0){
			this.nodeRows = subtreeReadIndexes;
		}
	}

	private void createChildren(){
		if(maxChildCount >= 1 && left == null){
			left = new TreeNode<T>(new Region(subtreeRows.start, nodeRows.start - 1), tree, this);
		}
		if(maxChildCount == 2 && right == null){
			right = new TreeNode<T>(new Region(nodeRows.end + 1, subtreeRows.end), tree, this);
		}
	}	

	private void fetchConcisedContent(AreaRequest areaRequest){
		tree.createFileRequest(areaRequest, this.nodeRows, this, areaRequest.status);
	}

	public void resetRequestDistributor(){
		requestDistributor = false;
		if(left != null){
			left.resetRequestDistributor();
		}
		if(right != null){
			right.resetRequestDistributor();
		}
	}


	public void processAreaRequest(AreaRequest areaRequest) {
		if(nodeBpRegion != null){

			//System.out.println(nodeBpRegion.getLength() + ", " +  areaRequest.getLength()  + ", " +  nodeRows.getLength());
			
			long nonZeroLength = Math.max(nodeBpRegion.getLength(), 1);
			
//			boolean detailNeeded = 
//				areaRequest.getLength() / nonZeroLength < RESOLUTION;
			
			boolean recursionNeeded = 
				areaRequest.getLength() / RESOLUTION < getSubtreeBpLengthApproximation();

			//Create file request or return value of this node
			if(nodeBpRegion.intercepts(areaRequest)){									
				
				if(areaRequest.status.concise){
					createConcisedResult(areaRequest.status);
				} else {
					//Concised value isn't enough, file has to be read
					
					tree.createFileRequest(areaRequest, nodeRows, this, areaRequest.status);
				}
			}				


			if(recursionNeeded){
				if(maxChildCount > 0 && left == null){
					createChildren();
				}

				if(left != null && areaRequest.start <= nodeBpRegion.start){ 
					left.processAreaRequest(areaRequest);
				}

				//Comparison to start, because datafile is sorted according to region start locations
				if(right != null && areaRequest.end >= nodeBpRegion.start){ 
					right.processAreaRequest(areaRequest);
				}	
			}

		} else {
			//Continue recursion when to data arrives
			//System.out.println("Area request delayed because of missing bp values " + areaRequest);
			fetchConcisedContent(areaRequest);
		}
	}
	
	private void createConcisedResult(FsfStatus status) {
		tree.createAreaResult(new AreaResult<RegionContent>(
				status, new RegionContent(nodeBpRegion, concisedValue)));
	}

	private long getSubtreeBpLengthApproximation(){
		if(parent != null && parent.nodeBpRegion != null && this.nodeBpRegion != null){
			return Math.abs(parent.nodeBpRegion.start - this.nodeBpRegion.start) * 2;
		} else {
			return Long.MAX_VALUE;
		}
	}

	/**
	 * @param fileResult null if should be calculated from children
	 */
	public void processFileResult(FileResult fileResult) {

		//if(fileResult != null){
		
		//System.out.println("ProcessFileResul");
		
		if(concisedValue == null || nodeBpRegion == null){
					
			ChunkParser parser = instructions.getParser();
			parser.setChunk(fileResult.chunk);
			
			concisedValue = instructions.conciser.concise(fileResult.chunk, parser);
			nodeBpRegion = new Region(
					parser.getLong(fileResult.chunk.readIndex, Content.BP_START),
					parser.getLong(
							(fileResult.chunk.readIndex + parser.getReadCount() - 1), 
							Content.BP_START));
													
			if(fileResult.request.areaRequest != null ){
//				System.out.println("Continue stopped area request: " + 
//						fileResult.request.areaRequest + 
//						", got node bp area: " + nodeBpRegion);
				
				//Continue stopped recursion now when bp location of this node is known				
				processAreaRequest(fileResult.request.areaRequest);
			}
		}
		
		if(nodeBpRegion.intercepts(fileResult.request.areaRequest)){
			//System.out.println(nodeBpRegion + ", " + concisedValue);
			if(fileResult.status.concise){
				createConcisedResult(fileResult.status);
			} else {				
				createReads(fileResult.chunk, fileResult.request.areaRequest, fileResult.status);				
			}
		}
	}
	
	public void createReads(ByteChunk chunk, AreaRequest areaRequest, FsfStatus status){
	//List<List<Object>> hits = new LinkedList<List<Object>>();
	ChunkParser parser = instructions.getParser();
	parser.setChunk(chunk);

	for(int i = 0; i < parser.getReadCount(); i++){	
		
		Region readReg = instructions.conciser.getRegion(parser, i + chunk.readIndex);
		
		
		Map<Content, Object> values = parser.getValues(
				i + chunk.readIndex, areaRequest.requestedContents);
		
		//if(areaRequest.requestedContents.contains(Content.FILE_INDEX)){
			values.put(Content.FILE_INDEX, i + chunk.readIndex);
		//}
		
				
		
//		
//		if("chr1".equals((String)read.get(instructions.fileDef.indexOf(Content.CHROMOSOME))) && 
//			region.intercepts(readReg)){
			
			//hits.add(read);
			tree.createAreaResult(new AreaResult<RegionContent>(
					status, new RegionContent(readReg, values)));
//		}
	}		
	//tree.createAreaResult(new AreaResult<List<List<Object>>>(status, hits, instructions.fileDef));
	}
}


