package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;

public class IndexKey  implements Comparable<IndexKey> {
	
	/**
	 * The natural order of these keys is primarily according to start positions. Field lineId 
	 * is the secondary sort condition so that lines with identical start position aren't lost and are 
	 * kept in original order.
	 * 
	 * @param start
	 * @param lineId Unique identifier for each line, preferably maintains the original order. For example
	 * line number or line start byte position.
	 */
	public IndexKey(BpCoord start, long lineId) {
		this.start = start;
		this.lineNumber = lineId;
	}
	
	private BpCoord start;
	private long lineNumber;

	@Override
	public boolean equals(Object o) {
		if (o instanceof IndexKey) {
			return (this.compareTo((IndexKey) o) == 0);
		}
		return false;
	}

	@Override
	public int hashCode() {
		return start.hashCode();
	}

	@Override
	public int compareTo(IndexKey o) {

		int startComparison = start.compareTo(o.start);

		if (startComparison != 0) {
			return startComparison;

		} else {
			return ((Long)lineNumber).compareTo(o.lineNumber);
		}
	}
}