package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

/**
 * Region of genome limited by two BpCoord locations. 
 *
 */
public class BpCoordRegion implements Comparable<BpCoordRegion> {
	public BpCoord start;
	public BpCoord end;

	public BpCoordRegion(BpCoord start, BpCoord end) {
		this.start = start;
		this.end = end;
	}

	public BpCoordRegion(Long start, Long end, Chromosome chr) {
		
		this.start = new BpCoord(start, chr);
		this.end = new BpCoord(end, chr);
	}

	public BpCoordRegion(Long start, Chromosome chr1, Long end, Chromosome chr2) {
		this.start = new BpCoord(start, chr1);
		this.end = new BpCoord(end, chr2);
	}

	public BpCoordRegion() {
		this(null, null);
	}

	public BpCoordRegion(BpCoordRegion bpRegion) {
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

	public BpCoordRegion clone() throws CloneNotSupportedException {
		return new BpCoordRegion(start, end);
	}

	public int compareTo(BpCoordRegion o) {
		
		int startComparison = start.compareTo(o.start);

		if (startComparison != 0) {
			return startComparison;
			
		} else {
			return end.compareTo(o.end);
		}
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof BpCoordRegion) {
			BpCoordRegion other = (BpCoordRegion) o;
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

	public boolean intersects(BpCoordRegion other) {
		BpCoord intersectionStart = start.max(other.start);
		BpCoord intersectionEnd = end.min(other.end);
		
		// Intersection has negative length <=> there is no intersection
		return intersectionStart.compareTo(intersectionEnd) <= 0;
	}

	public BpCoordRegion intersect(BpCoordRegion other) {
		if (!intersects(other)) {
			throw new IllegalArgumentException("regions do not intersect");
		}
		return new BpCoordRegion(start.max(other.start), end.min(other.end));
	}

	public BpCoordRegion merge(BpCoordRegion other) {
		if (!intersects(other)) {
			throw new IllegalArgumentException("regions do not intersect");
		}
		return new BpCoordRegion(start.min(other.start), end.max(other.end));
	}
}
