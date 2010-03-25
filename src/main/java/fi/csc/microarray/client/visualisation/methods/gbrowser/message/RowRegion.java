package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

public class RowRegion implements Comparable<RowRegion> {
	public Long start;
	public Long end;

	public RowRegion(Long start, Long end) {
		this.start = start;
		this.end = end;
	}

	public RowRegion() {
		this(null, null);
	}

	public double getLength() {
		return end - start;
	}

	public double getMid() {
		return (start + end) / 2.0;
	}

	public String toString() {
		return "RowRegion [" + start + " - " + end + "]";
	}

	public RowRegion clone() {
		return new RowRegion(start, end);
	}

	public boolean intercepts(RowRegion other) {
		return end >= other.start && start < other.end;
	}

	public RowRegion intercept(RowRegion other) {
		return new RowRegion(Math.max(start, other.start), Math.min(end, other.end));
	}

	public int compareTo(RowRegion o) {
		int startComparison = start.compareTo(start);

		if (startComparison != 0) {
			return startComparison;
			
		} else {
			return end.compareTo(o.end);
		}
	}

	@Override
	public boolean equals(Object o) {
		
		if (o instanceof RowRegion) {
			RowRegion other = (RowRegion) o;
			return start.equals(other.start) && end.equals(other.end);
		}
		return false;
	}

	@Override
	public int hashCode() {
		return start.hashCode();
	}

	public boolean contains(Long point) {
		return point >= start && point < end;
	}
}
