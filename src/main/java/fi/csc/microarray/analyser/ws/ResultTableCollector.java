package fi.csc.microarray.analyser.ws;

import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.LinkedList;
import java.util.Set;

public class ResultTableCollector {
	
	public static class ResultRow {
		
		private HashMap<String, ResultField> fields = new HashMap<String, ResultField>();

		public void addValue(String fieldName, String value) {
			if (fields.containsKey(fieldName)) {
				throw new IllegalStateException(fieldName + " already added");
			}
			fields.put(fieldName, new ResultField(value));
		}
		
		public String getValue(String fieldName) {
			if (!fields.containsKey(fieldName)) {
				throw new IllegalArgumentException("field " + fieldName + " not found");			
			}
			return fields.get(fieldName).getValue();
		}
	}
	
	private static class ResultField {
		private String value;
		
		public ResultField(String value) {
			this.value = value;
		}

		public String getValue() {
			return value;
		}
	}
	
	private LinkedList<ResultRow> resultRows = new LinkedList<ResultRow>();
	private Set<String> fieldNames = new LinkedHashSet<String>();
	
	
	public void addAnnotation(int index, String fieldName, String value) {
		addRowsTo(index);
		resultRows.get(index).addValue(fieldName, value);
		fieldNames.add(fieldName);
	}
	
	
	public static interface RowFilter {
		public boolean shouldRemove(ResultRow row);
	}
	
	public void filterRows(RowFilter filter) {
		LinkedList<ResultRow> removed = new LinkedList<ResultRow>();
		for (ResultRow row : resultRows) {
			if (filter.shouldRemove(row)) {
				removed.add(row);
			}
		}
		resultRows.removeAll(removed);
	}

	
	public String[][] asTable(String[] columns) {
		String[][] table = new String[resultRows.size()][columns.length];
		int n = 0;
		int m = 0;
		for (ResultRow row : resultRows) {
			for (String column : columns) {
				table[n][m] = row.getValue(column);
				m++;
			}
			n++;
			m = 0;
		}
		return table;
	}	

	private void addRowsTo(int index) {
		while (resultRows.size() < (index+1)) {
			resultRows.add(new ResultRow());
		}
	}
	
	public String[] getFieldNames() {
		return fieldNames.toArray(new String[fieldNames.size()]);
	}

}
