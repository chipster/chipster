package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class Gene extends TreeSet<RegionContent> implements Comparable<Gene>{
	
	public BpCoordRegion region;
	public String id;
	
	public Gene(BpCoordRegion region, String id) {
		this.region = region;
		this.id = id;
	}

	@Override
	public boolean add(RegionContent part){
		
		return super.add(part);
	}

	@Override
	public int hashCode(){
		return id.hashCode();
	}

	@Override
	public boolean equals(Object o){
		Gene other = (Gene) o;		
		return id.equals(other.id);
	}

	public int compareTo(Gene other) {
		return region.compareTo(other.region);
	}
}
