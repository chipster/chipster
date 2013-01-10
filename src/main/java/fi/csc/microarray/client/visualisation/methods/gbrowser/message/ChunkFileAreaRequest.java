package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.Collection;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;

/**
 * A request for content from the given view area. The view layer uses these requests to get data from the processing layer.  
 *
 */
public class ChunkFileAreaRequest extends AreaRequest {

	private static int MAX_RECURSION_DEPTH = 8;
	
	public int depthToGo;
	
	private ChunkFileStatus status;
	
	/* 
	 * Return ChunkFileStatus instead of @link DataRetrievalStatus  
	 */
	@Override
	public DataRetrievalStatus getStatus() {
		return status;
	}

	/**
	 * Constructs a new request with max depth available (no splits done yet).
	 *  
	 * @see #AreaRequest(Region, Collection, DataRetrievalStatus, int)
	 */
	public ChunkFileAreaRequest(Region region, Collection<ColumnType> requestedContents, DataRetrievalStatus status) {
		this(region, requestedContents, status, MAX_RECURSION_DEPTH);
	}
	
	/**
	 * Constructs a new request. 
	 * 
	 * @param region 				coordinate region for which content is needed
	 * @param requestedContents		what fields of the content file are needed
	 * @param status				status				
	 * @param depthToGo				counter for counting tree split depth
	 */
	public ChunkFileAreaRequest(Region region, Collection<ColumnType> requestedContents, DataRetrievalStatus status, int depthToGo) {
		super(region, requestedContents, status);
		this.depthToGo = depthToGo;
	}
	
	/**
	 * Clone clones also relevant fields so that clones are independent.
	 */
	@Override
	public ChunkFileAreaRequest clone() throws CloneNotSupportedException {
		return new ChunkFileAreaRequest(this, this.getRequestedContents(), this.getStatus().clone(), this.depthToGo);
	}
}
