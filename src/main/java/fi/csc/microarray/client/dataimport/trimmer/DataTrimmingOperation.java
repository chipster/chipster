package fi.csc.microarray.client.dataimport.trimmer;

/**
 * Class that represents a single data modification operatio.
 * 
 * @author mkoski
 *
 */
public abstract class DataTrimmingOperation {

	public static final int ALL_COLUMNS = Integer.MAX_VALUE;
	private int columnIndex;
	
	/**
	 * Creates a new DataTrimmingOperation
	 * 
	 * @param oldString Old string
	 * @param newString New string
	 * @param columnIndex index of the column to which the operation is affected
	 * @param isRegexp is the operation done by using regular expressions
	 */
	public DataTrimmingOperation(int columnIndex) {
		this.columnIndex = columnIndex;
	}

	public int getColumnIndex() {
		return columnIndex;
	}
	
	public abstract String doTrimming(String stringToTrim);
	
}
