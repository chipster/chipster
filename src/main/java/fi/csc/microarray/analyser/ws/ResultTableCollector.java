package fi.csc.microarray.analyser.ws;

import java.util.HashMap;
import java.util.LinkedList;

public class ResultTableCollector {
	
	public static class ResultRow {
		
		private HashMap<String, ResultField> fields = new HashMap<String, ResultField>();

		public void addValue(String fieldName, String value) {
			if (fields.containsKey(fieldName)) {
				throw new IllegalStateException(fieldName + " already added");
			}
			fields.put(fieldName, new ResultField(fieldName, value));
		}
		
		public String getValue(String fieldName) {
			if (!fields.containsKey(fieldName)) {
				throw new IllegalArgumentException("field " + fieldName + " not found");			
			}
			return fields.get(fieldName).getValue();
		}
	}
	
	private static class ResultField {
		private String name;
		private String value;
		
		public ResultField(String name, String value) {
			this.name = name;
			this.value = value;
		}

		public String getName() {
			return name;
		}

		public String getValue() {
			return value;
		}

		public void setName(String name) {
			this.name = name;
		}

		public void setValue(String value) {
			this.value = value;
		}
		
	}
	
	private LinkedList<ResultRow> resultRows = new LinkedList<ResultRow>();
	
	public void addAnnotation(int index, String fieldName, String value) {
		addRowsTo(index);
		resultRows.get(index).addValue(fieldName, value);
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

}
