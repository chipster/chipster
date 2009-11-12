package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.Region;

public class ConstantLengthChunker extends Chunker {
	
	private ConstantLengthChunkParser parser;
	
	public ConstantLengthChunker(ConstantLengthChunkParser parser){
		this.parser = parser;
	}

	@Override
	public int getChunkByteLength() {
		return (int)parser.getRowLength() * 32;
	}

	@Override
	public long getFilePosition(long readIndex) {
		return readIndex * parser.getRowLength();
	}

	@Override
	public long getReadIndex(long filePosition) {
		return filePosition / parser.getRowLength();
	}
	
	@Override
	public Region getChunkMiddleOf(Region readIndexes){
		
		Region nodeRows = new Region();
		
		//Round to next chunk split
		nodeRows.start = (long)Math.ceil(readIndexes.getMid() / getChunkReadLength()) * getChunkReadLength();
		nodeRows.end = nodeRows.start + getChunkReadLength() - 1;
		
		return nodeRows;
	}
	
	@Override
	public int getChildCount(Region subtreeReadIndexes){
		
		if(subtreeReadIndexes.getLength() <= getChunkReadLength()){
			return 0;
		} else if(subtreeReadIndexes.getLength() <= getChunkReadLength() * 2){
			return 1;
		} else {
			return 2;
		}
	}

	@Override
	public int getChunkReadLength() {
		return (int) (getChunkByteLength() / parser.getRowLength());
	}
}
