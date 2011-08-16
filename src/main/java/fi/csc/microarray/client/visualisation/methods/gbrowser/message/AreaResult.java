package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

/**
 * Result with content for some view area. The tree layers uses these results to send content back to view layer.
 *
 */
public class AreaResult {

	public FsfStatus status;
	public RegionContent content;

	public AreaResult(FsfStatus status, RegionContent content) {
		this.status = status;
		this.content = content;
	}
}
