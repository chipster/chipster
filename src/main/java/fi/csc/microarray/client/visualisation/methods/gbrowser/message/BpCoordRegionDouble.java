package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

/**
 * Floating-point version of BpCoordRegion. Used when floating point precision is needed in calculations.
 * 
 * @see BpCoordRegion
 *
 */
public class BpCoordRegionDouble implements Comparable<BpCoordRegionDouble> {
	
	public BpCoordDouble start;
	public BpCoordDouble end;

	public BpCoordRegionDouble(BpCoordDouble start, BpCoordDouble end) {
		this.start = start;
		this.end = end;
	}

	public BpCoordRegionDouble(Double start, Double end, Chromosome chr) {
		this.start = new BpCoordDouble(start, chr);
		this.end = new BpCoordDouble(end, chr);
	}

	public BpCoordRegionDouble(Double start, Chromosome chr1, Double end, Chromosome chr2) {
		this.start = new BpCoordDouble(start, chr1);
		this.end = new BpCoordDouble(end, chr2);
	}

	public BpCoordRegionDouble() {
		this(null, null);
	}

	public BpCoordRegionDouble(BpCoordRegion reg) {
		this((double)reg.start.bp, reg.start.chr, (double)reg.end.bp, reg.end.chr);
	}

	public double getLength() {
		return end.minus(start);
	}

	public double getMid() {
		return start.plus(end) / 2;
	}

	public String toString() {
		return "Region [" + start + " - " + end + "]";
	}

	public BpCoordRegionDouble clone() {
		return new BpCoordRegionDouble(start, end);
	}

	public boolean intercepts(BpCoordRegionDouble other) {
		return end.compareTo(start) > 0 && start.compareTo(end) < 0;
	}

	public BpCoordRegionDouble intercept(BpCoordRegionDouble other) {
		return new BpCoordRegionDouble(start.max(other.start), end.min(other.end));
	}

	public int compareTo(BpCoordRegionDouble o) {
		int startComparison = start.compareTo(start);

		if (startComparison != 0) {
			return startComparison;
		} else {
			return end.compareTo(o.end);
		}
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof BpCoordRegionDouble) {
			BpCoordRegionDouble other = (BpCoordRegionDouble) o;
			return start.equals(other.start) && end.equals(other.end);
		}
		return false;
	}

	@Override
	public int hashCode() {
		return start.hashCode();
	}

	public boolean contains(BpCoordDouble point) {
		return point.compareTo(start) >= 0 && point.compareTo(end) < 0;
	}

	public void move(double bpMove) {
		start.bp += bpMove;
		end.bp += bpMove;
	}

}
