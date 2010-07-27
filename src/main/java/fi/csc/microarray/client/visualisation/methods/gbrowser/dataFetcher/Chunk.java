package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

public class Chunk implements Cloneable {
	
	private long byteLocation;
	private String content;
	
	public Chunk() {		
	}

	public Chunk(String content) {
		this.content = content;
	}

	public Chunk(String content2, long byteLocation2) {
		this.content = content2;
		this.byteLocation = byteLocation2;
	}

	public String getContent() {
		return content;
	}

	public void setContent(String content) {
		this.content = content;
	}

	public long getByteLocation() {
		return byteLocation;
	}

	public void setByteLocation(long byteLocation) {
		this.byteLocation = byteLocation;
	} 
	
	public Chunk clone() {
		return new Chunk(new String(content), byteLocation);
	}
}
