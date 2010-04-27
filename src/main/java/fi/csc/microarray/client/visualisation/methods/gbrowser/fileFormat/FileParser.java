package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Collection;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
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
	 * @return position in actual file, that is contentBytePosition offset by header
	 */
	public long getFilePosition(long contentBytePosition) {
		return contentBytePosition;
	}
}
