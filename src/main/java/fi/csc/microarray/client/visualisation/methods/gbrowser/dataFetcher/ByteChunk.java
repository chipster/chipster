package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

@Deprecated
public class ByteChunk {

	public byte[] byteContent;

	public ByteChunk(int length) {

		this.byteContent = new byte[length];
	}

}