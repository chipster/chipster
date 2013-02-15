package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;


/**
 * Result with content for some view area. The processing layer uses these results to send content back to view layer.
 *
 */
public class AreaResult {

	private DataRetrievalStatus status;
	private List<RegionContent> contents;

	public AreaResult(DataRetrievalStatus status, List<RegionContent> contents) {
		this.status = status;
		this.contents = contents;
	}
	
	public AreaResult(DataRetrievalStatus status, Region region, Object value) {
		this.status = status;
		
		LinkedHashMap<ColumnType, Object> valueMap = new LinkedHashMap<ColumnType, Object>();
		
		valueMap.put(ColumnType.VALUE, value);
		
		RegionContent regionContent = new RegionContent(region, valueMap);
		
		List<RegionContent> list = new LinkedList<RegionContent>();
		list.add(regionContent);
		
		this.contents = list;
	}

	public DataRetrievalStatus getStatus() {
		return status;
	}

	public List<RegionContent> getContents() {
		return contents;
	}

}
