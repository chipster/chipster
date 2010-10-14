package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

public class CigarItem {
	private Long length;
	private String type;
	
	public CigarItem(long length, String type) {
		this.length = length;
		this.type = type;
	}
	
	public long getLength() {
		return length;
	}
	
	public String getType() {
		return type;
	}
}
