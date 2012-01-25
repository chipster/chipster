package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.LinkedHashMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Type;

/**
 * Content for given genomic region. Content is data dependent, but basically it is data parsed from the 
 * input file or other data source. All the rows fall within the genomic region.
 *
 */
public class RegionContent implements Comparable<RegionContent> {
	
	public BpCoordRegion region;
	public LinkedHashMap<ColumnType, Object> values;

	public RegionContent(BpCoordRegion region, LinkedHashMap<ColumnType, Object> values) {
		this.region = region;
		this.values = values;
	}

	public RegionContent(BpCoordRegion region, Object concisedValue) {
		this.region = region;
		this.values = new LinkedHashMap<ColumnType, Object>();
		this.values.put(ColumnType.VALUE, concisedValue);
	}

	public RegionContent(BpCoordRegion region, Object concisedValueForward,  Object concisedValueReverse) {
		this.region = region;
		this.values = new LinkedHashMap<ColumnType, Object>();
		this.values.put(ColumnType.VALUE_FORWARD, concisedValueForward);
		this.values.put(ColumnType.VALUE_REVERSE, concisedValueReverse);
	}

	public RegionContent(BpCoordRegion region) {
		this.region = region;
		this.values = new LinkedHashMap<ColumnType, Object>(); 
	}

	public int compareTo(RegionContent other) {

		int regionComparison = this.region.compareTo(other.region);

		if (regionComparison != 0) {
			return regionComparison;			
			
		} else {
			return values.toString().compareTo(other.values.toString());
		}
	}
	
	@Override
	public int hashCode() {
		return region.hashCode();
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof RegionContent) {
			return this.compareTo((RegionContent) o) == 0;
		} else {
			return false;
		}
	}
	
	@Override
	public String toString() {
		String extra = "";
		for (Object value : values.values()) {
			extra += "\t" + Type.toString(value);
		}
		return region.toString(true) + extra;
	}
}
