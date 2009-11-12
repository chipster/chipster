package fi.csc.microarray.client.visualisation.methods.genomeBrowser.message;

import java.util.HashMap;
import java.util.Map;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.Content;


public class RegionContent implements Comparable<RegionContent> {
	public Region region;
	public Map<Content, Object> values;
	
	public RegionContent(Region region, Map<Content, Object> values){
		this.region = region;
		this.values = values;
	}

	public RegionContent(Region region, Object concisedValue) {
		this.region = region;
		this.values = new HashMap<Content, Object>();
		this.values.put(Content.VALUE, concisedValue);
	}

	public int compareTo(RegionContent other) {
		
		int regionComparison = this.region.compareTo(other.region);
		if(regionComparison != 0){
			return regionComparison;
		}
		
		Long first = (Long)values.get(Content.FILE_INDEX);
		Long second = (Long)other.values.get(Content.FILE_INDEX);
		
		if(first == null && second == null){
			return 0;
		}

		if(second == null){
			return 1;
		}
		
		if(first == null){
			return -1;
		}
		
		return first.compareTo(second);	
	}
	
//	@Override
//	public int hashCode(){
//		return region.hashCode();
//	}
//	
//	@Override
//	public boolean equals(Object o){
//		RegionContent other = (RegionContent) o;		
//		return region.equals(other.region);
//	}
}
