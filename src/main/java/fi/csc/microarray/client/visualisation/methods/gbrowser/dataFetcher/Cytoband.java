package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;

public class Cytoband implements Comparable<Cytoband> {
	
	private BpCoordRegion region;
	private String band;
	private Stain stain; 
	
	public enum Stain {

		GNEG("gneg"), 
		GPOS25("gpos25"),
		GPOS33("gpos33"),
		GPOS50("gpos50"),
		GPOS66("gpos66"),
		GPOS75("gpos75"), 
		GPOS100("gpos100"), 
		ACEN("acen"), 
		GVAR("gvar"),
		STALK("stalk"),
		TIP("tip"),
		UNRECOGNIZED(null);

		private String id;

		Stain(String id) {
			this.id = id;
		}

		public String getId() {
			return id;
		}
	}

	public Cytoband(BpCoordRegion region, String band, String stain) {
		this.region = region;
		this.band = band;
		this.stain = getStain(stain);
	}
	
	public Cytoband(BpCoordRegion region, String band, Stain stain) {
		this.region = region;
		this.band = band;
		this.stain = stain;
	}
	
	public Cytoband(BpCoord coord) {
		this.region = new BpCoordRegion(coord, coord);
	}

	private Stain getStain(String id) {
		for (Stain stain : Stain.values()) {
			if (stain.getId().equals(id)) {
				return stain;
			}
		}
		return Stain.UNRECOGNIZED;
	}
	
	

	public int compareTo(Cytoband other) {

		int comparison = this.region.compareTo(other.region);

		if (comparison != 0) {
			return comparison;			
		} 
		
		comparison = this.band.compareTo(other.getBand());
		
		if (comparison != 0 && this.band != null) {
			return comparison;			
		} 
		
		if (this.stain != null) {
			comparison = this.stain.compareTo(other.getStain());
		}
		
		return comparison;
	}
	
	public Stain getStain() {
		return stain;
	}

	public String getBand() {
		return band;
	}

	@Override
	public int hashCode() {
		return region.hashCode();
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof Cytoband) {
			return this.compareTo((Cytoband) o) == 0;
		} else {
			return false;
		}
	}
	
	@Override
	public String toString() {
		
		return region.toString(true) + ", " + band + ", " + stain.getId();
	}

	public BpCoordRegion getRegion() {
		return region;
	}
}
