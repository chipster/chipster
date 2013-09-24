package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.text.DecimalFormat;
import java.util.LinkedHashMap;


/**
 * Content for given genomic region. Content is data dependent, but basically it is data parsed from the 
 * input file or other data source. All the rows fall within the genomic region.
 *
 */
public class Feature implements Comparable<Feature> {
	
	public Region region;
	public LinkedHashMap<DataType, Object> values;

	public Feature(Region region, LinkedHashMap<DataType, Object> values) {
		this.region = region;
		this.values = values;
	}

	public Feature(Region region) {
		this.region = region;
		this.values = new LinkedHashMap<DataType, Object>(); 
	}

	public int compareTo(Feature other) {

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
		if (o instanceof Feature) {
			return this.compareTo((Feature) o) == 0;
		} else {
			return false;
		}
	}
	
	@Override
	public String toString() {
		
		final DecimalFormat FLOAT_FORMAT = new DecimalFormat("0.#########");
		String extra = "";
		
		for (Object value : values.values()) {
			
			if (value instanceof Float) {
				extra += "\t" + FLOAT_FORMAT.format((Float)value);
				
			} else {
				extra += "\t" + value.toString();
			}
		}
		return region.toString(true) + extra;
	}

	public IndexKey getIndexKey() {
		return (IndexKey)values.get(DataType.ID);
	}

	public Object getValueObject() {
		return values.get(DataType.VALUE);
	}
}
