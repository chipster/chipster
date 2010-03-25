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
		if (regionComparison != 0) {
			return regionComparison;
		}

		Long first = (Long) values.get(ColumnType.FILE_INDEX);
		Long second = (Long) other.values.get(ColumnType.FILE_INDEX);

		if (first == null && second == null) {
			return 0;
		}

		if (second == null) {
			return 1;
		}

		if (first == null) {
			return -1;
		}

		return first.compareTo(second);
	}
}
