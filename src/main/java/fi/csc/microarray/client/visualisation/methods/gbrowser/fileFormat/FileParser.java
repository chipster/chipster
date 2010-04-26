package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Collection;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ByteRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public abstract class FileParser {

	public String chunk;

	public abstract List<RegionContent> getAll(Collection<ColumnType> requestedContents);	
	public abstract RegionContent[] concise(BpCoordRegion readIndexRegion);
	/**
	 * @return Region where start is the start position of the first row of chunk and the end is 
	 * the start position of the last row of chunk 
	 */
	public abstract BpCoordRegion getBpRegion();
	public abstract String getName();
	public abstract FileParser clone();
	public abstract long getDefaulChunkLength();

	public void setChunk(String chunk) {
		this.chunk = chunk;
	}
	
	public void clone(FileParser copy) {
		copy.chunk = new String(this.chunk);
	}
	
	/**
	 * @param position in content part of the file
	 * @return position in actual file
	 */
	public long getFilePosition(long contentBytePosition) {
		return contentBytePosition;
	}

	public ByteRegion getChunkRegionMiddleOf(ByteRegion byteRegion) {

		ByteRegion nodeRegion = new ByteRegion();

		// round to next chunk split
		nodeRegion.start = (long) Math.ceil(byteRegion.getMid() / getDefaulChunkLength()) * getDefaulChunkLength();
		
		Long defaultEnd = nodeRegion.start + getDefaulChunkLength() - 1;
		nodeRegion.end = Math.min(defaultEnd, byteRegion.end);

		return nodeRegion;
	}

	public int getChildCount(ByteRegion subtreeByteRegion) {

		if (subtreeByteRegion.getLength() <= getDefaulChunkLength()) {
			return 0;
		} else if (subtreeByteRegion.getLength() <= getDefaulChunkLength() * 2) {
			return 1;
		} else {
			return 2;
		}
	}
}
