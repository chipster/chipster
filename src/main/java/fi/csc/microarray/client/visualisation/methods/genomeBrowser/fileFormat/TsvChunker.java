package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.Region;

public class TsvChunker extends Chunker {

	@Override
	public int getChunkByteLength() {
		return 1024*8;
	}

	@Override
	public long getFilePosition(long readIndex) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public long getReadIndex(long filePosition) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public int getChildCount(Region subtreeReadIndexes) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public Region getChunkMiddleOf(Region subtreeReadIndexes) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public int getChunkReadLength() {
		// TODO Auto-generated method stub
		return 0;
	}
}
