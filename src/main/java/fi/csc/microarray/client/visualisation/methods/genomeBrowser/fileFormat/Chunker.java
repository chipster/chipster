package fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.Region;

public abstract class Chunker {

	public abstract long getFilePosition(long readIndex);
	
	public abstract long getReadIndex(long filePosition);
	
	public abstract int getChunkByteLength();
	
	public abstract int getChunkReadLength();

	public abstract Region getChunkMiddleOf(Region subtreeReadIndexes);

	public abstract int getChildCount(Region subtreeReadIndexes);
}
