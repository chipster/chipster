package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.HashMap;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;

/**
 * Content for given genomic region. Content is data dependent, but basically it is data parsed from tabular data. All the 
 * rows fall within the genomic region.
 *
 */
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
		RegionContent other = (RegionContent) o; 
		return region.equals(other.region) && values.toString().equals(other.values.toString());		
	}
}
