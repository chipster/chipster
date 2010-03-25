package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

public class ByteChunk {

	public long rowIndex;
	public long byteLength;
	public byte[] byteContent;

	public ByteChunk(int length) {
		this.rowIndex = -1;
		this.byteLength = -1;
		this.byteContent = new byte[length];
	}

}