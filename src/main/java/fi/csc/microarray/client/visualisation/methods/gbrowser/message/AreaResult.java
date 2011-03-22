package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

/**
 * Result with content for some view area. The tree layers uses these results to send content back to view layer.
 *
 */
public class AreaResult<T> {

	public FsfStatus status;
	public T content;

	public AreaResult(FsfStatus status, T content) {
		this.status = status;
		this.content = content;
	}
}
