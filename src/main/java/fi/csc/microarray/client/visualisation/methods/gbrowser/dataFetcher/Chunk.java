package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

/**
 * DOCME
 * 
 * Chunk is a part of file containing only full lines. The default length of the chunk can be 
 * defined in the parser classes, but the exact length will be only in the some order of 
 * magnitude and the exact length depends on how the lines happen to become divided between 
 * these chunks. 
 * 
 * @author klemela
 */
public class Chunk implements Cloneable {

	// Location of the first byte of chunk in the original data file
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
