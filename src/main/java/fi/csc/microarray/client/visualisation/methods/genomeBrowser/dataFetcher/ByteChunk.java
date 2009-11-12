package fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher;
	public class ByteChunk{
		
		public ByteChunk(int length){
			this.readIndex = -1;
			this.length = -1;
			content = new byte[length];;
		}
		
		public long readIndex;
		public long length;
		public byte[] content;
	}