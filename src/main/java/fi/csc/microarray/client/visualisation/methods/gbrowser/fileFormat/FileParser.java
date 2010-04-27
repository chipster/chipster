package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Collection;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public abstract class FileParser {

	public abstract List<RegionContent> getAll(String chunk, Collection<ColumnType> requestedContents);	
	public abstract RegionContent[] concise(String chunk);
	/**
	 * @return Region where start is the start position of the first row of chunk and the end is 
	 * the start position of the last row of chunk 
	 */
	public abstract BpCoordRegion getBpRegion(String chunk);
	public abstract String getName();
	public abstract long getDefaulChunkLength();
	
	/**
	 * @param position in content part of the file
	 * @return position in actual file, that is contentBytePosition offset by header
	 */
	public long getFilePosition(long contentBytePosition) {
		return contentBytePosition;
	}
}
