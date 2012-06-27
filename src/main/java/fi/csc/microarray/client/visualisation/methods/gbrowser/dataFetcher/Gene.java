package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.Collection;
import java.util.HashMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

public class Gene extends HashMap<String, Transcript> implements Comparable<Gene> {
	
	private Region region;
	private String name;
	private String biotype;
	private String id;

	public Gene(String name, String biotype, String id) {
		this.name = name;
		this.biotype = biotype;
		this.id = id;
	}

	public int compareTo(Gene other) {
		return this.id.compareTo(other.id);
	}

	@SuppressWarnings("unused") // might be needed in future
	private String getBiotype() {
		return biotype;
	}

	@Override
	public int hashCode() {
		return id.hashCode();
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

	public void addTranscript(String transcId, Transcript transcript) {
		
		this.put(transcId, transcript);
	}


	public Collection<Transcript> getTranscripts() {
		return this.values() ;
	}

	public String getName() {
		return name;
	}

	public void addExon(Exon exon, String geneId, String transcId, String transcName) {
		Transcript transc;
		
		if ((transc = this.get(transcId)) == null) {
			transc = new Transcript(transcName, this, transcId);
			this.put(transcId, transc);
		}
		
		transc.addExon(exon, transcId);
		exon.setTranscript(transc);
		
		if(region == null) {
			try {
				region = transc.getRegion().clone();
			} catch (CloneNotSupportedException e) {
				e.printStackTrace(); //Shouldn't happen
			}
		} else {
			
			this.region = region.fill(transc.getRegion());
		}
	}

	public void setRegion(Region region) {
		this.region = region;
	}
}
