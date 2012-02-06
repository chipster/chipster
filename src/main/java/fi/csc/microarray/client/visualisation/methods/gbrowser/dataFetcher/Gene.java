package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.LinkedList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

public class Gene implements Comparable<Gene> {
	
	private Region region;
	private List<Transcript> tempTranscripts = new LinkedList<Transcript>();
	private SortedSet<Transcript> transcripts;
	private String name;
	private String biotype;

	public Gene(String name, String biotype) {
		this.name = name;
		this.biotype = biotype;
	}

	public int compareTo(Gene other) {

		int regionComparison = this.region.compareTo(other.region);
		int biotypeComparison = 0;
		
		if (regionComparison != 0) {
			return regionComparison;
		}
		
		if (this.biotype != null) {
			biotypeComparison = this.biotype.compareTo(other.getBiotype());
		} else if (other.getBiotype() != null) {
			biotypeComparison = 1;
		}
	
		return biotypeComparison;
	}

	private String getBiotype() {
		return biotype;
	}

	@Override
	public int hashCode() {
		return region.hashCode();
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof Gene) {
			return this.compareTo((Gene) o) == 0;
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

	public void addTranscript(Transcript transcript) {
		
		tempTranscripts.add(transcript);
		
		if(region == null) {
			try {
				region = transcript.getRegion().clone();
			} catch (CloneNotSupportedException e) {
				e.printStackTrace(); //Shouldn't happen
			}
		} else {
			
			this.region = region.fill(transcript.getRegion());
		}
	}

	public void prepareForReading() {
		
		transcripts = new TreeSet<Transcript>(tempTranscripts);
		tempTranscripts = null;
	}

	public SortedSet<Transcript> getTranscripts() {
		return transcripts;
	}

	public String getName() {
		return name;
	}
}
