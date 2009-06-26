package fi.csc.microarray.databeans.features;

import java.util.LinkedHashMap;
import java.util.Map;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.databeans.DataBean;

public class ConstantTableFeature extends BasicFeature {

	public class ConstantFloatTable implements Table {

		private int row = -1;
		private Map<String, Integer> columnMap;
		private Object[][] values;
		
		public ConstantFloatTable(Map<String, Integer> columnMap, Object[][] values) {
			this.columnMap = columnMap;
			this.values = values;
		}

		public boolean nextRow() {
			row++;
			return row < values[0].length;
		}

		public String[] getColumnNames() {
			return columnMap.keySet().toArray(new String[0]);
		}

		public float getFloatValue(String columnName) {
			return (Float)getValue(columnName);
		}

		public String getStringValue(String columnName) {
			return (String)getValue(columnName);
		}

		public int getIntValue(String columnName) {
			return (Integer)getValue(columnName);
		}

		public Object getValue(String columnName) {
			return values[columnMap.get(columnName)][row];
		}

		public boolean hasColumn(String columnName) {
			return columnMap.get(columnName) != null;
		}

		public int getColumnCount() {
			return columnMap.keySet().size();
		}

		public void close() {
			// we don't need to do anything			
		}
	}

	private Map<String, Integer> columnMap = new LinkedHashMap<String, Integer>();
	private Object[][] values;
	
	public ConstantTableFeature(DataBean bean, FeatureProvider factory, String[] columns, Object[][] values) {
		super(bean, factory);
		if (columns.length != values.length) {
			throw new IllegalArgumentException("columns and values do not match in size");
		}
		for (int i = 0; i < columns.length; i++) {
			columnMap.put(columns[i], i);
		}
		this.values = values;
	}

	public Table asTable() throws MicroarrayException {
		return new ConstantFloatTable(columnMap, values);
	}
}
