package fi.csc.microarray.databeans.features.table;

import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;

import org.apache.log4j.Logger;

public class EditableTable {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(EditableTable.class);
	
	private final static String SEPARATOR = "\t";
	
	private HashMap<String, List<String>> columns = new HashMap<String, List<String>>();
	/**
	 * For indexing into actual columns.
	 */
	private LinkedList<String> columnNames = new LinkedList<String>();
	
	public void addColumns(Map<String, List<String>> newColumns) {
		columnNames.addAll(newColumns.keySet());
		columns.putAll(newColumns);
	}
	
	/**
	 * Adds a new columns after the specified column index.
	 * 
	 * @see #addColumn(String, List)
	 */
	public void addColumn(String columnName, int afterColumn, List<String> values) {
		
		if (columns.containsKey(columnName)) {
			throw new IllegalArgumentException("column " + columnName + " exists already");
		}
		
		columnNames.add(afterColumn, columnName);
		columns.put(columnName, values);
	}

	public void removeColumn(String columnName) {
		columnNames.remove(columnName);
		columns.remove(columnName);
	}
	
	public void addColumn(String columnName, List<String> values) {
		
		if (columns.containsKey(columnName)) {
			throw new IllegalArgumentException("column " + columnName + " exists already");
		}
		
		columnNames.add(columnName);
		columns.put(columnName, values);
	}

	public int getColumnCount() {
		return columns.size();
	}

	public void setValue(String columnName, int rowIndex, String value) {
		columns.get(columnName).set(rowIndex, value);
	}
	
	/**
	 * Return row count. Checks that all columns have same row count.
	 * 
	 * @return number of rows in the matrix
	 * @throws RuntimeException if counts do not match
	 */
	public int getRowCount() {
		int count = -1;
		for (List<String> column : columns.values()) {
			if (count == -1) {
				count = column.size();
			}
			if (count != column.size()) {
				throw new RuntimeException("column counts " + count + " and " + column.size() + " do not match");
			}
			
		}
		return count;		
	}

	public void writeTo(OutputStream out) {
		logger.debug("writing data out, has " + columns.keySet().size() + " columns");
		
		// write column names
		PrintWriter writer = new PrintWriter(new OutputStreamWriter(out)); 
		Iterator<String> columnIter = columnNames.iterator();
		while (columnIter.hasNext()) {
			String name = columnIter.next();
			writer.print(name);
			if (columnIter.hasNext()) {
				writer.print(SEPARATOR);
			}
		}
		writer.print("\n");
		
		// write data
		for (int i = 0; i < getRowCount(); i++) {
			Iterator<String> columnIter2 = columnNames.iterator();
			while (columnIter2.hasNext()) {
				List<String> column = columns.get(columnIter2.next());
				writer.print(column.get(i));
				if (columnIter2.hasNext()) {
					writer.print(SEPARATOR);
				}
			}
			writer.print("\n");
		}
		writer.flush();
		logger.debug("" + getRowCount() + " rows written");
	}

	public Iterable<String> getColumnNames() {
		return columnNames;
	}

	public String getValue(String columnName, int rowIndex) {
		return columns.get(columnName).get(rowIndex);
	}

	public String getColumnName(int columnIndex) {
		return columnNames.get(columnIndex);		
	}

	public boolean containsColumn(String columnName) {
		return columnNames.contains(columnName);
	}
}
