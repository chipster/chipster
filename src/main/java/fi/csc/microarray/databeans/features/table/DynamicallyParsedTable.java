package fi.csc.microarray.databeans.features.table;

import java.io.BufferedReader;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.LinkedList;

import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.databeans.features.table.TableColumnProvider.MatrixParseSettings;
import fi.csc.microarray.databeans.features.table.TableColumnProvider.TableColumn;
import fi.csc.microarray.util.LookaheadLineReader;


/**
 * Implements all actual tabular data parsing.  
 * 
 * @author Aleksi Kallio
 *
 */
public class DynamicallyParsedTable implements Table {

	private LookaheadLineReader source;
	private boolean headerParsed;
	private MatrixParseSettings settings;
	private LinkedList<Integer> columnNumbers;
	private HashMap<String, String> values;
	private String[] columnNames;
	private DataBean bean;

	public DynamicallyParsedTable(DataBean bean, MatrixParseSettings settings, LinkedList<Integer> columnNumbers) {
		this.bean = bean;
		this.settings = settings;
		this.columnNumbers = columnNumbers;
		this.columnNames = settings.columns.keySet().toArray(new String[0]);
		reset();
	}

	/**
	 * Reset table iteration. After call to reset iteration will start from the first row.
	 */
	public void reset() {
		try {
			this.source = new LookaheadLineReader(new BufferedReader(new InputStreamReader(bean.getContentByteStream())));
			this.headerParsed = false;

		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}
	
	public boolean nextRow() {
		try {

			// we stop at 1) EOS (null), 2) footer starter and 3) empty row (if header is parsed, because header can contains empty rows)
			if (source.peekLine() == null || 
					(settings.footerStarter != null && source.peekLine().contains(settings.footerStarter)) ||
					(headerParsed && "".equals(source.peekLine().trim()))) {
				
				values = null; // trying to read values will result now in error
				return false; // signal that we are at end
			}

			// do we have to first process header?
			if (!headerParsed) {
				
				// parse away headers, if any
				if (settings.headerTerminator != null) {
					TableColumn.parseAwayHeader(source, settings);
				}

				// skip column name row, if any
				if (settings.hasColumnNames) {
					source.readLine();
				}
				headerParsed = true;
			}
			
			if (source.peekLine() == null) {
				return false; // header parsing has eaten all content
			}
			
			ArrayList<String> row = parseRow(source.readLine());
			values = new HashMap<String, String>();
			
			if (columnNumbers.size() == 0) {
				// return all
				int i = 0; 
				for (String string : row) {
					values.put(columnNames[i++], string);
				}
				
			} else {
				for (Integer number : columnNumbers) {
					values.put(columnNames[number], row.get(number));
				}
			}

			return true; // not at end yet
			
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}

	private String preprocessExternalData(String string) {
		return string.replace("EMPTY", "NaN");
	}
	
	private ArrayList<String> parseRow(String row) throws MicroarrayException {

		ArrayList<String> result = new ArrayList<String>(settings.columns.size());
		row = preprocessExternalData(row);

		String[] cells = TableColumn.tokeniseRow(row);
		for (int i = 0; i < settings.columns.values().size(); i++) {
			String cell;
			if (i < cells.length) {
				cell = cells[i];
			} else {
				// we are stuffing too short rows with empty cells (should this worry us?)
				cell = "";
			}
			result.add(cell);
		}
		
		return result;
	}

	public String[] getColumnNames() {
		String[] colunmNameSlice = new String[columnNumbers.size()];
		for (int i = 0; i < columnNumbers.size(); i++) {
			colunmNameSlice[i] = columnNames[columnNumbers.get(i)];
		}
		return colunmNameSlice;
	}

	public float getFloatValue(String columnName) {
		try {
			return new Float(values.get(columnName));
		} catch (NumberFormatException nfe) {
			return new Float(Float.NaN);
		} catch (NullPointerException ne) {
			throw new IllegalArgumentException("column name " + columnName + " was not found");
		}
	}

	public int getIntValue(String columnName) {
		return (int)getFloatValue(columnName);
	}

	public String getStringValue(String columnName) {
		return values.get(columnName);
	}

	public Object getValue(String columnName) {
		try {
			return new Float(values.get(columnName));
		} catch (NumberFormatException e) {
			return getStringValue(columnName);
		}
	}

	public boolean hasColumn(String columnName) {
		for (String name : columnNames) {
			if (name.equals(columnName)) {
				return true;
			}
		}
		return false;
	}

	public int getColumnCount() {
		return columnNames.length;
	}
}
