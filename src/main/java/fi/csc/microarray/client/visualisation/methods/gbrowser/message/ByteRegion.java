package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

public class ByteRegion implements Comparable<ByteRegion> {
	public Long start;
	public Long end;
	
	
	/**
	 * If the particular file position isn't read yet, the exact location of the row can't be known.
	 * False value in this attribute means that the region doesn't start or end in the positions of 
	 * start and end attributes, but just after the next new line characters following both locations.
	 * There is one exception to this  interpretation: start position of 0 must not be moved to
	 * following new line character, otherwise the first line of file will be lost.  
	 */
	public boolean exact;

	public ByteRegion(Long start, Long end, boolean exact) {
		this.start = start;
		this.end = end;
		this.exact = exact; 
	}

	public ByteRegion() {
		this(null, null, false);
	}

	public long getLength() {
		return end - start;
	}

	public double getMid() {
		return (start + end) / 2.0;
	}

	public String toString() {
		return "RowRegion [" + start + " - " + end + "]";
	}

	public ByteRegion clone() {
		return new ByteRegion(start, end, exact);
	}

	public boolean intercepts(ByteRegion other) {
		return end >= other.start && start < other.end;
	}

	public ByteRegion intercept(ByteRegion other) {
		return new ByteRegion(Math.max(start, other.start), Math.min(end, other.end), exact);
	}

	public int compareTo(ByteRegion o) {
		int startComparison = start.compareTo(start);

		if (startComparison != 0) {
			return startComparison;
			
		} else {
			return end.compareTo(o.end);
		}
	}

	@Override
	public boolean equals(Object o) {
		
		if (o instanceof ByteRegion) {
			ByteRegion other = (ByteRegion) o;
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
