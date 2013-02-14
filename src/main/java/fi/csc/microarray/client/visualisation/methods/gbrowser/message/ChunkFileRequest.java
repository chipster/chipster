package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.TreeNode;

/**
 * Request for file content. Used by the processing (tree) layer to request content from the data source layer. 
 *
 */
public class ChunkFileRequest {

	public ChunkFileAreaRequest areaRequest;
	public ByteRegion byteRegion;
	public TreeNode node;

	public ChunkFileStatus status;

	public ChunkFileRequest(ChunkFileAreaRequest areaRequest, ByteRegion rowRegion, TreeNode node, ChunkFileStatus status) {
		super();
		this.byteRegion = rowRegion;
		this.node = node;
		this.status = status;
		this.areaRequest = areaRequest;
	}
}