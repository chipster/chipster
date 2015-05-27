package fi.csc.microarray.databeans.features;

/**
 * Interface for table like iterable constructs with named columns. 
 * Quite similar to JDBC ResultSet.
 * 
 * @author akallio
 *
 */
public interface Table extends AutoCloseable {

	/**
	 * Iterates to next table row and return true iff we are still on a valid
	 * row.
	 */
	public boolean nextRow();
	
	/**
	 * Returns Float representation of the value or NaN if conversion fails.
	 */
	public float getFloatValue(String columnName);
	
	/**
	 * Return String representation of the value.
	 */
	public String getStringValue(String columnName);

	/**
	 * Returns a int representation of the value or NaN if conversion fails.
	 */
	public int getIntValue(String columnName);

	/**
	 * Returns a Float representation if value is numeric, String otherwise.
	 */
	public Object getValue(String columnName);
	
	public String[] getColumnNames();

	public boolean hasColumn(String columnName);

	/**
	 * Returns number of columns.
	 */
	public int getColumnCount();
	
	/**
	 * Closes underlying resources and makes the table not usable any more.
	 */
	public void close();
}
