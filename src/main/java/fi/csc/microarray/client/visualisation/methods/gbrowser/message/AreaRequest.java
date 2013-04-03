package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.Collection;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;

/**
 * A request for content from the given view area. The view layer uses these requests to get data from the processing layer.  
 *
 */
public class AreaRequest extends Region {
	
	private DataRetrievalStatus status;
	private Collection<ColumnType> requestedContents;
	
	/**
	 * Constructs a new request. 
	 * 
	 * @param region 				coordinate region for which content is needed
	 * @param requestedContents		what fields of the content file are needed
	 * @param status				status				
	 * @param depthToGo				counter for counting tree split depth
	 */
	public AreaRequest(Region region, Collection<ColumnType> requestedContents, DataRetrievalStatus status) {
		super(region.start, region.end);
		this.requestedContents = requestedContents;
		this.status = status;
	}
	
	/**
	 * Clone clones also relevant fields so that clones are independent.
	 */
	@Override
	public AreaRequest clone() throws CloneNotSupportedException {
		return new AreaRequest(this, this.requestedContents, this.getStatus().clone());
	}
	
	public DataRetrievalStatus getStatus() {
		return status;
	}
	
	public Collection<ColumnType> getRequestedContents() {
		return requestedContents;
	}
}
