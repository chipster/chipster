package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

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

	public String toString() {
		return "Region [" + start + " - " + end + "]";
	}

	public BpCoordRegion clone() throws CloneNotSupportedException {
		return new BpCoordRegion(start, end);
	}

	public boolean intercepts(BpCoordRegion other) {
		return other.end.compareTo(start) > 0 && other.start.compareTo(end) < 0;
	}

	public BpCoordRegion intercept(BpCoordRegion other) {
		return new BpCoordRegion(start.max(other.start), end.min(other.end));
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
		return point.compareTo(start) >= 0 && point.compareTo(end) < 0;
	}
}
