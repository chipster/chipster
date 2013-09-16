package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.Collection;
import java.util.TreeSet;


public class Transcript extends TreeSet<Exon> implements Comparable<Transcript> {
	
	private Gene gene;
	private Region region;
	private String name;
	private String id;

	public Transcript(String name, Gene gene, String id) {
		this.name = name;
		this.gene = gene;
		this.id = id;
	}

	public int compareTo(Transcript other) {

		
		return this.id.compareTo(other.id);
//		int regionComparison = this.region.compareTo(other.region);
//		
//		return regionComparison;
	}

	@Override
	public int hashCode() {
		return id.hashCode();
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof Transcript) {
			return this.compareTo((Transcript) o) == 0;
		} else {
			return false;
		}
	}
	
	@Override
	public String toString() {
		
		return region.toString(true) + ", " + name;
	}

	public Region getRegion() {
		return region;
	}

	public void addExon(Exon exon, String transcId) {
		
		this.add(exon);
		
		if(region == null) {
			try {
				region = exon.getRegion().clone();
			} catch (CloneNotSupportedException e) {
				e.printStackTrace();
			}
		} else {
			
			this.region = region.fill(exon.getRegion());
		}
	}

	public String getName() {
		return name;
	}
	
	public String getId() {
		return id;
	}
	public Collection<Exon> getExons() {
		return this;
	}
}
