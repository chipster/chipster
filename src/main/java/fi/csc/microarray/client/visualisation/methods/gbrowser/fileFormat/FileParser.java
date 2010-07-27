package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.Collection;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Generic file parser. Information is divided into chunks.
 * 
 * @See FsfStatus for description of chunk
 */
public abstract class FileParser {


	/**
	 * DOCME
	 * 
	 * Method for getting all the lines from the data chunk. RequestedContents is used to define
	 * which columns should be returned to avoid parsing data that isn't really used. It's up to
	 * the parser implementation to decide how these ColumnTypes are mapped to the columns of 
	 * particular file format. Of course the Track that will eventually draw to content must be 
	 * able to support the data type that is returned.  
	 * 
	 * @See FsfStatus for description of chunk
	 * @param chunk
	 * @param requestedContents
	 * @return list of RegionContent objects containing the requested columns from the chunk data
	 */
	public abstract List<RegionContent> getAll(String chunk, Collection<ColumnType> requestedContents);
	
	/**
	 * Method for getting a summary of data in the chunk. These summaries will be used, when the 
	 * View area is too big to be able to draw all the detailed content in that area. Good summary
	 * is for example read coverage, which is: 
	 * 
	 * count of reads in chunk * read length / length of total chunk area
	 * 
	 * This can be also be weighted with the read quality to calculate quality coverage.
	 * 
	 * Even if the parser doesn't support concised results, it must still return an empty array.
	 * 
	 * @See FsfStatus for description of chunk
	 * @param chunk
	 * @return 
	 */
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
