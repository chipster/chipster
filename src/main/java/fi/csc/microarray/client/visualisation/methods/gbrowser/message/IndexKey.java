package fi.csc.microarray.client.visualisation.methods.gbrowser.message;


/**
 * IndexKey object is an identifier of the data row in data file. A same line of the data
 * produces always equal IndexKey objects and different lines produce always differing IndexKeys. The natural order 
 * of the IndexKey objects is the same as the original order of the lines in the file.
 * 
 * Optionally the IndexKey may contain to row number of the line.
 * 
 * @author klemela
 */
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
		this.lineId = lineId;
	}
	
	/**
	 * Same as the contructor with two parameters, but sets the flag isRowNumber. When this flag is true, the 
	 * lineId must be the line number of the line in file. In other words, lineId of the first line is 0, lineId 
	 * of the second line is 1 and so on.  
	 * 
	 * @param start
	 * @param lineId
	 * @param isRowIndex
	 */
	public IndexKey(BpCoord start, long lineId, boolean isRowNumber) {
		this(start, lineId);
		this.isRowNumber = isRowNumber;
	}
	
	private BpCoord start;
	private long lineId;
	private boolean isRowNumber = false;

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
			return ((Long)lineId).compareTo(o.lineId);
		}
	}
	
	@Override
	public String toString() {
		return "IndexKey [Start: " + start + " Id: " + lineId + "]";
	}

	/**
	 * Get the row number of the line in the original data file, excluding
	 * all comment lines. The row nubmer of first content line is 0.
	 * 
	 * @return
	 */
	public Integer getRowNumber() {
		if (isRowNumber) {
			return (int) lineId;
		}
		return null;
	}
}