package fi.csc.microarray.client.visualisation.methods.gbrowser.message;


/**
 * A class with similar features with {@link IndexKey} offering an identifier and natural order for the data items.
 * While IndexKey is an identifier of the original data line, this class can be used with higher level objects
 * where the use of data line identifiers is not sensible. For example, this is the case when the item is constructed from 
 * multiple lines.   
 * 
 * @author klemela
 */
public class PositionAndStringKey  implements Comparable<PositionAndStringKey> {
	
	/**
	 * The natural order of these keys is primarily according to start positions. Field lineId 
	 * is the secondary sort condition so that lines with identical start position aren't lost and are 
	 * kept in original order.
	 * 
	 * @param start
	 * @param lineId Unique identifier for each line, preferably maintains the original order. For example
	 * line number or line start byte position.
	 */
	public PositionAndStringKey(BpCoord start, String lineId) {
		this.start = start;
		this.lineId = lineId;
	}
	
	private BpCoord start;
	private String lineId;

	@Override
	public boolean equals(Object o) {
		if (o instanceof PositionAndStringKey) {
			return (this.compareTo((PositionAndStringKey) o) == 0);
		}
		return false;
	}

	@Override
	public int hashCode() {
		return start.hashCode();
	}

	@Override
	public int compareTo(PositionAndStringKey o) {

		int startComparison = start.compareTo(o.start);

		if (startComparison != 0) {
			return startComparison;

		} else {
			return ((String)lineId).compareTo(o.lineId);
		}
	}
	
	@Override
	public String toString() {
		return "IndexKey [Start: " + start + " Id: " + lineId + "]";
	}
}