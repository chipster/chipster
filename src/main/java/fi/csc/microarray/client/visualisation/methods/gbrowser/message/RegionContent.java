package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.HashMap;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;

public class RegionContent implements Comparable<RegionContent> {
	
	public BpCoordRegion region;
	public Map<ColumnType, Object> values;

	public RegionContent(BpCoordRegion region, Map<ColumnType, Object> values) {
		this.region = region;
		this.values = values;
	}

	public RegionContent(BpCoordRegion region, Object concisedValue) {
		this.region = region;
		this.values = new HashMap<ColumnType, Object>();
		this.values.put(ColumnType.VALUE, concisedValue);
	}

	public int compareTo(RegionContent other) {

		int regionComparison = this.region.compareTo(other.region);

		return regionComparison;
	}
}
