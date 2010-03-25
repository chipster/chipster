package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

public class BpCoordDouble implements Comparable<BpCoordDouble> {

	public Double bp;
	public Chromosome chr;

	public BpCoordDouble(Double bp, Chromosome chr) {
		this.bp = bp;
		this.chr = chr;
	}

	public BpCoordDouble(BpCoord coord) {
		this.bp = (double) (long) coord.bp;
		this.chr = coord.chr;
	}

	public int compareTo(BpCoordDouble o) {

		int chrComparison = chr.compareTo(((BpCoordDouble) o).chr);

		if (chrComparison != 0) {
			return chrComparison;
		} else {
			return bp.compareTo(((BpCoordDouble) o).bp);
		}
	}

	public boolean equals(Object o) {
		if (o instanceof BpCoordDouble) {
			BpCoordDouble other = (BpCoordDouble) o;

			return chr.equals(other.chr) && bp.equals(other.bp);
		}
		return false;
	}

	public int hashCode() {
		return bp.hashCode();
	}

	public double minus(BpCoordDouble o) {
		if (chr.equals(o.chr)) {
			return bp - o.bp;
		}
		return Double.NaN;
	}

	public double plus(BpCoordDouble o) {
		if (chr.equals(o.chr)) {
			return bp + o.bp;
		}
		return Double.NaN;
	}

	public BpCoordDouble move(double increment) {
		return new BpCoordDouble(bp + increment, chr);
	}

	public String toString() {
		return "Bp: " + bp + ", Chr: " + chr;
	}

	public BpCoordDouble max(BpCoordDouble other) {
		return this.compareTo(other) < 0 ? other : this;
	}

	public BpCoordDouble min(BpCoordDouble other) {
		return this.compareTo(other) > 0 ? other : this;
	}

	public BpCoord asBpCoord() {
		return new BpCoord((long) (double) bp, chr);
	}

}