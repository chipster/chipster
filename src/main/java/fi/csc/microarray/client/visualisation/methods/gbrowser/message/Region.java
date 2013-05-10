package fi.csc.microarray.client.visualisation.methods.gbrowser.message;


/**
 * Region of genome limited by two {@link BpCoord} coordinates. 
 *
 */
public class Region implements Comparable<Region> {
	public BpCoord start;
	public BpCoord end;
	
	@Deprecated
	public Strand strand;

	public Region(BpCoord start, BpCoord end) {
		this.start = start;
		this.end = end;
	}
	
	@Deprecated
	public Region(BpCoord start, BpCoord end, Strand strand) {
		this.start = start;
		this.end = end;
		this.strand = strand;
	}
	
	@Deprecated
	public Region(Long start, Long end, Chromosome chr, Strand strand) {
		this.start = new BpCoord(start, chr);
		this.end = new BpCoord(end, chr);
		this.strand = strand;
	}

	public Region(Long start, Long end, Chromosome chr) {
		
		this.start = new BpCoord(start, chr);
		this.end = new BpCoord(end, chr);
	}

	public Region(Long start, Chromosome chr1, Long end, Chromosome chr2) {
		this.start = new BpCoord(start, chr1);
		this.end = new BpCoord(end, chr2);
	}

	public Region() {
		this(null, null);
	}

	public Region(Region bpRegion) {
		this(new BpCoord(bpRegion.start), new BpCoord(bpRegion.end));
	}

	public Long getLength() {
		return end.minus(start);
	}

	public Long getMid() {
		return start.plus(end) / 2;
	}

	@Override
	public String toString() {
		return toString(false);
	}
	
	public String toString(boolean tsvFormat) {
		if (tsvFormat) {
			return start.chr + "\t" + start.bp + "\t" + end.bp;
		} else {
			return "Region [" + start + " - " + end + "]";
		}
	}

	public Region clone() throws CloneNotSupportedException {
		return new Region(start.clone(), end.clone(), strand);
	}

	public int compareTo(Region o) {
		
		int startComparison = start.compareTo(o.start);

		if (startComparison != 0) {
			return startComparison;
			
		} else {
			return end.compareTo(o.end);
		}
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof Region) {
			Region other = (Region) o;
			return start.equals(other.start) && end.equals(other.end);
		}
		return false;
	}

	@Override
	public int hashCode() {
		return start.hashCode();
	}

	public boolean contains(BpCoord point) {
		return start.chr.equals(point.chr) && point.compareTo(start) >= 0 && point.compareTo(end) < 0;
	}
	
	public boolean contains(Long point) {
		return point.compareTo(start.bp) >= 0 && point.compareTo(end.bp) < 0;
	}

	/**
	 * Return true if this region intersects with the other. Handles chromosomes as sequential. i.e. if the region ends are in different chromosomes, is 
	 * region interpreted to include:
	 *  - everything after the start position in the chromosome of start position
	 *  - everything before the end position in the chromosome of end position
	 *  - all chromosomes that are between these chromosomes in the order defined by Chromosome.compareTo method
	 * 
	 * @param other
	 * @return
	 */
	public boolean intersects(Region other) {
	
		BpCoord intersectionStart = start.max(other.start);
		BpCoord intersectionEnd = end.min(other.end);
		
		// Intersection has negative length <=> there is no intersection
		return intersectionStart.compareTo(intersectionEnd) <= 0;
	}

	public Region intersect(Region other) {
		if (!intersects(other)) {
			throw new IllegalArgumentException("regions do not intersect");
		}
		return new Region(start.max(other.start), end.min(other.end));
	}

	public Region merge(Region other) {
		if (!intersects(other)) {
			throw new IllegalArgumentException("regions do not intersect");
		}
		return new Region(start.min(other.start), end.max(other.end));
	}
	
	public Region fill(Region other) {
		
		BpCoord left = start.min(end.min(other.start.min(other.end)));
		BpCoord right = start.max(end.max(other.start.max(other.end)));

		return new Region(left, right, this.strand);
	}

	public Strand getStrand() {
		return strand;
	}

	public Region grow(long growLength) {
		
		double center = (start.bp + end.bp) / 2;
		long length = getLength() + growLength;
	
		return new Region((long) (center - length / 2), (long) (center + length / 2), start.chr);  
	}
}
