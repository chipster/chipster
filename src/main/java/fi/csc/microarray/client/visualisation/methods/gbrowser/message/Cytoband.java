package fi.csc.microarray.client.visualisation.methods.gbrowser.message;


public class Cytoband implements Comparable<Cytoband> {
	
	private Region region;
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
		GPOS("gpos"),
		ACEN("acen"), 
		GVAR("gvar"),
		STALK("stalk"),
		TIP("tip"),
		WHITE("white"),
		UNRECOGNIZED(null);

		private String id;

		Stain(String id) {
			this.id = id;
		}

		public String getId() {
			return id;
		}
	}

	public Cytoband(Region region, String band, String stain) {
		this.region = region;
		this.band = band;
		this.stain = getStain(stain);
	}
	
	public Cytoband(Region region, String band, Stain stain) {
		this.region = region;
		this.band = band;
		this.stain = stain;
	}
	
	public Cytoband(BpCoord coord) {
		this.region = new Region(coord, coord);
	}

	private Stain getStain(String id) {
		for (Stain stain : Stain.values()) {
			if(stain.getId() != null && stain.getId().equals(id)) {
				return stain;
			}
		}
		return Stain.UNRECOGNIZED;
	}
	
	

	public int compareTo(Cytoband other) {
		
		int regionComparison = this.region.compareTo(other.region);
		
		int bandComparison = 0;
		int stainComparison = 0;
		
		if (this.band != null && other.getBand() != null) {
			bandComparison = this.band.compareTo(other.getBand());
		} else if (this.band == null && other.getBand() == null) {
			bandComparison = 0;
		} else if (this.band == null) {
			bandComparison = -1;
		}
		
		if (this.stain != null && other.getBand() != null) {
			stainComparison = this.stain.compareTo(other.getStain());
		} else if (this.stain == null && other.getBand() == null) {
			stainComparison = 0;
		} else if (this.stain == null) {
			stainComparison = -1;
		}
		
		if (regionComparison != 0) {
			return regionComparison;
		} else if (bandComparison != 0) {
			return bandComparison;
		}
		return stainComparison;
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

	public Region getRegion() {
		return region;
	}
}
