package fi.csc.microarray.client.visualisation.methods.genomeBrowser.track;

import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.RegionContent;

public class Gene extends TreeSet<RegionContent> implements Comparable<Gene>{	

	public BpCoord minBp = BpCoord.getMax();
	public BpCoord maxBp = BpCoord.getMin();
	public String id = null;	

	@Override
	public boolean add(RegionContent part){

		minBp = minBp.min(part.region.start);
		maxBp = maxBp.max(part.region.end);
		id  = (String)part.values.get(ColumnType.ID);
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
