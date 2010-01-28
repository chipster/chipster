package fi.csc.microarray.databeans.features;

import java.util.List;

import fi.csc.microarray.exception.MicroarrayException;

/**
 * @author akallio
 */
public class RestrictModifier implements Modifier {

	public static final int RESTRICT_TO_ROWS = 2000;
	
	/**
	 * Delegates to original table, but restricts row iteration.
	 *
	 */
	private static class RestrictedTable implements Table {

		private Table table;
		private int row = 0;
		
		public RestrictedTable(Table table) {
			this.table = table;
		}
		
		public String[] getColumnNames() {
			return table.getColumnNames();
		}

		public int getColumnCount() {
			return table.getColumnCount();
		}

		public float getFloatValue(String columnName) {
			return table.getFloatValue(columnName);
		}

		public int getIntValue(String columnName) {
			return table.getIntValue(columnName);
		}

		public String getStringValue(String columnName) {
			return table.getStringValue(columnName);
		}

		public Object getValue(String columnName) {
			return table.getValue(columnName);
		}

		public boolean hasColumn(String columnName) {
			return table.hasColumn(columnName);
		}

		public boolean nextRow() {
			this.row++;
			return row < RESTRICT_TO_ROWS && table.nextRow();
		}

		public void close() {
			table.close();			
		}
		
	}
	private static class RestrictModifierFeature extends ModifiedFeature {

		private Feature original;

		protected RestrictModifierFeature(List<Feature> inputs) {
			super(inputs);
			if (inputs.size() != 1) {
				throw new IllegalArgumentException("restrict modifier must have 1 parameter");
			}
			this.original = inputs.get(0);
		}

		public Table asTable() throws MicroarrayException {
			Table table = original.asTable();
			if (table == null) {
				return null;
			} else {
				return new RestrictedTable(table);
			}
		}
	}

	private RestrictModifierFeature output;
	
	public Feature getOutput() {
		return output;
	}

	public void setInputs(List<Feature> inputs) {
		this.output = new RestrictModifierFeature(inputs);
	}

}
