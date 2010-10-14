package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.TreeNode;

/**
 * Request for file content. Used by the tree layer to request content from the file layer. 
 *
 */
public class FileRequest {

	public AreaRequest areaRequest;
	public ByteRegion byteRegion;
	public TreeNode node;

	public FsfStatus status;

	public FileRequest(AreaRequest areaRequest, ByteRegion rowRegion, TreeNode node, FsfStatus status) {
		super();
		this.byteRegion = rowRegion;
		this.node = node;
		this.status = status;
		this.areaRequest = areaRequest;
	}

}