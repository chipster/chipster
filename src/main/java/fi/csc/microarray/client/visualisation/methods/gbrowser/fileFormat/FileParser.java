package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Collection;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Generic file parser. Information is divided into chunks.
 * 
 * DOCME: what is a chunk?
 *
 */
public abstract class FileParser {

    /**
     * DOCME
     */
	public abstract List<RegionContent> getAll(String chunk, Collection<ColumnType> requestedContents);
	
	public abstract RegionContent[] concise(String chunk);
	
	/**
	 * @return Region where start is the start position of the first row of chunk and the end is 
	 * the start position of the last row of chunk 
	 */
	public abstract BpCoordRegion getBpRegion(String chunk);

	/**
	 * @return Human readable name of the file format.
	 */
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
