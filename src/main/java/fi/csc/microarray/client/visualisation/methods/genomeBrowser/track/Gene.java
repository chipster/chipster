package fi.csc.microarray.client.visualisation.methods.genomeBrowser.track;

import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.Content;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.RegionContent;

public class Gene extends TreeSet<RegionContent> implements Comparable<Gene>{	

	public Long minBp = Long.MAX_VALUE;
	public long maxBp = 0;
	public String id = null;	

	@Override
	public boolean add(RegionContent part){

		minBp = Math.min(minBp, part.region.start);
		maxBp = Math.max(maxBp, part.region.end);
		id  = (String)part.values.get(Content.ID);
		return super.add(part);
	}

	@Override
	public int hashCode(){
		return minBp.hashCode();
	}

	@Override
	public boolean equals(Object o){
		Gene other = (Gene) o;		
		return minBp.equals(other.minBp);
	}

	public int compareTo(Gene other) {
		return minBp.compareTo(other.minBp);
	}

}
