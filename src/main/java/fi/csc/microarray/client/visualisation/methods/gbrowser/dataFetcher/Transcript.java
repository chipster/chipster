package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.SortedSet;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

public class Transcript implements Comparable<Transcript> {
	
	private Region region;
	private SortedSet<Exon> exons = new TreeSet<Exon>();
	private String name;

	public Transcript(String name) {
		this.name = name;
	}

	public int compareTo(Transcript other) {

		int regionComparison = this.region.compareTo(other.region);
		
		return regionComparison;
	}

	@Override
	public int hashCode() {
		return region.hashCode();
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

	public void addExon(Exon exon) {
		exons.add(exon);
		
		if(region == null) {
			try {
				region = exon.getRegion().clone();
			} catch (CloneNotSupportedException e) {
				e.printStackTrace(); //Shouldn't happen
			}
		} else {
			
			this.region = region.fill(exon.getRegion());
		}
	}

	public String getName() {
		return name;
	}

	public SortedSet<Exon> getExons() {
		return exons;
	}
}
