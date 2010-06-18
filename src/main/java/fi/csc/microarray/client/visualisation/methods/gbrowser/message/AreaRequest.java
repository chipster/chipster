package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.Collection;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;

/**
 * A request for content from the given view area. The view layers uses these requests to get data from the tree layer.  
 *
 */
public class AreaRequest extends BpCoordRegion {

	private static int MAX_RECURSION_DEPTH = 8;
	
	public FsfStatus status;
	public Collection<ColumnType> requestedContents;
	public int depthToGo;

	/**
	 * Constructs a new request with max depth available (no splits done yet).
	 *  
	 * @see #AreaRequest(BpCoordRegion, Collection, FsfStatus, int)
	 */
	public AreaRequest(BpCoordRegion region, Collection<ColumnType> requestedContents, FsfStatus status) {
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
	public AreaRequest(BpCoordRegion region, Collection<ColumnType> requestedContents, FsfStatus status, int depthToGo) {
		super(region.start, region.end);
		this.requestedContents = requestedContents;
		this.status = status;
		this.depthToGo = depthToGo;
	}
	
	/**
	 * Clone clones also relevant fields so that clones are independent.
	 */
	@Override
	public AreaRequest clone() throws CloneNotSupportedException {
		return new AreaRequest(this, this.requestedContents, this.status.clone(), this.depthToGo);
	}
}
