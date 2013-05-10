package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.LinkedList;
import java.util.List;



/**
 * Result with content for some view area. The processing layer uses these results to send content back to view layer.
 *
 */
public class DataResult {

	private DataStatus status;	
	private List<RegionContent> contents;
	private DataRequest request;

	public DataResult(DataStatus status, List<RegionContent> contents) {
		this.status = status;
		this.contents = contents;
	}	

	public DataResult(DataRequest request,
			LinkedList<RegionContent> contents) {
		this(request.getStatus(), contents);
		this.request = request;
	}

	public DataStatus getStatus() {
		return status;
	}

	public List<RegionContent> getContents() {
		return contents;
	}

	public DataRequest getRequest() {
		return request;
	}

}
