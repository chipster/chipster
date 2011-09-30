package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

/**
 * Basepair coordinate: base(pair) location and chromosome.
 *
 */
public class BpCoord implements Comparable<BpCoord> {

	private static final BpCoord MAX = new BpCoord(Long.MAX_VALUE, null);
	private static final BpCoord MIN = new BpCoord(0l, null);
	public Long bp;
	public Chromosome chr;

	public BpCoord(Long bp, Chromosome chr) {
		this.bp = bp;
		this.chr = chr;
	}

	public BpCoord(BpCoord o) {
		this(new Long(o.bp), new Chromosome(o.chr));
	}

	public BpCoord clone() {
		return new BpCoord(bp, chr);
	}

	public int compareTo(BpCoord o) {
		
		int chrComparison = chr.compareTo(o.chr);

		if (chrComparison != 0) {
			return chrComparison;
		} else {
			return bp.compareTo(o.bp);
		}
	}
	
	public int compareTo(Long o) {

		return bp.compareTo(o);		
	}

	public boolean equals(Object o) {
		if (o instanceof BpCoord) {
			BpCoord other = (BpCoord) o;

			return chr.equals(other.chr) && bp.equals(other.bp);
		}
		return false;
	}

	public int hashCode() {
		return bp.hashCode();
	}

	public Long minus(BpCoord o) {
		
		if (chr.equals(o.chr)) {
			return bp - o.bp;
		}
		return null;
	}

	public Long plus(BpCoord o) {
		if (chr.equals(o.chr)) {
			return bp + o.bp;
		}
		return null;
	}

	public String toString() {
		return "Bp: " + bp + ", Chr: " + chr;
	}

	public BpCoord max(BpCoord other) {
		return this.compareTo(other) < 0 ? other : this;
	}

	public BpCoord min(BpCoord other) {
		return this.compareTo(other) > 0 ? other : this;
	}

	public static BpCoord getMax() {
		return MAX.clone();
	}

	public static BpCoord getMin() {
		return MIN.clone();
	}
}