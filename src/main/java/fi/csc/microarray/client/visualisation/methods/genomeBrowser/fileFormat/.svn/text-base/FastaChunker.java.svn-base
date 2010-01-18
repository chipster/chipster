package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;


public class FastaChunker extends ConstantLengthChunker {

	private FastaChunkParser parser;

	public FastaChunker(FastaChunkParser parser) {
		super(null);
		
		this.parser = parser;
	}

	@Override
	public int getChunkByteLength() {
		return (int) parser.getRowLength() * 128;
	}

	@Override
	public long getFilePosition(long readIndex) {
		return readIndex * parser.getRowLength() + parser.getTitleLength();
	}

	@Override
	public long getReadIndex(long filePosition) {
		return (filePosition - parser.getTitleLength()) / parser.getRowLength();
	}

	@Override
	public int getChunkReadLength() {
		return (int) getReadIndex(getChunkByteLength());
	}
}
