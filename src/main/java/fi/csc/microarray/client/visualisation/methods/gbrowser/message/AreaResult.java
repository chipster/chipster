package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.List;

/**
 * Result with content for some view area. The processing layer uses these results to send content back to view layer.
 *
 */
public class AreaResult {

	private FsfStatus status;
	private List<RegionContent> contents;

	public AreaResult(FsfStatus status, List<RegionContent> contents) {
		this.status = status;
		this.contents = contents;
	}

	public FsfStatus getStatus() {
		return status;
	}

	public List<RegionContent> getContents() {
		return contents;
	}

}
