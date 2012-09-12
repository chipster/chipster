package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

/**
 * Floating-point version of {@link Region}. Used when floating point precision is needed in calculations.
 * 
 */
public class RegionDouble implements Comparable<RegionDouble> {
	
	public BpCoordDouble start;
	public BpCoordDouble end;

	public RegionDouble(BpCoordDouble start, BpCoordDouble end) {
		this.start = start;
		this.end = end;
	}

	public RegionDouble(Double start, Double end, Chromosome chr) {
		this.start = new BpCoordDouble(start, chr);
		this.end = new BpCoordDouble(end, chr);
	}

	public RegionDouble(Double start, Chromosome chr1, Double end, Chromosome chr2) {
		this.start = new BpCoordDouble(start, chr1);
		this.end = new BpCoordDouble(end, chr2);
	}

	public RegionDouble() {
		this(null, null);
	}

	public RegionDouble(Region reg) {
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

	public RegionDouble clone() {
		return new RegionDouble(start, end);
	}

	public boolean intercepts(RegionDouble other) {
		return end.compareTo(start) > 0 && start.compareTo(end) < 0;
	}

	public RegionDouble intercept(RegionDouble other) {
		return new RegionDouble(start.max(other.start), end.min(other.end));
	}

	public int compareTo(RegionDouble o) {
		int startComparison = start.compareTo(start);

		if (startComparison != 0) {
			return startComparison;
		} else {
			return end.compareTo(o.end);
		}
	}

	@Override
	public boolean equals(Object o) {
		if (o instanceof RegionDouble) {
			RegionDouble other = (RegionDouble) o;
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
