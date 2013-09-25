package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.Collection;


/**
 * A request for content from the given view area. The view layer uses these requests to get data from the processing layer.  
 *
 */
public class DataRequest extends Region {
	
	private DataStatus status;
	private Collection<DataType> requestedContents;
	
	/**
	 * Constructs a new request. 
	 * 
	 * @param region 				coordinate region for which content is needed
	 * @param requestedContents		what fields of the content file are needed
	 * @param status				status				
	 * @param depthToGo				counter for counting tree split depth
	 */
	public DataRequest(Region region, Collection<DataType> requestedContents, DataStatus status) {
		super(region.start, region.end);
		this.requestedContents = requestedContents;
		this.status = status;
	}
	
	/**
	 * Clone clones also relevant fields so that clones are independent.
	 */
	@Override
	public DataRequest clone() throws CloneNotSupportedException {
		return new DataRequest(this, this.requestedContents, new DataStatus(this.getStatus()));
	}
	
	public DataStatus getStatus() {
		return status;
	}
	
	public Collection<DataType> getRequestedContents() {
		return requestedContents;
	}
}
