package fi.csc.microarray.client.visualisation.methods.genomeBrowser.message;


public class AreaResult<T> {
	
	public FsfStatus status;
	public T content;
	
	public AreaResult(FsfStatus status, T content) {
		this.status = status;
		this.content = content;
	}
}
