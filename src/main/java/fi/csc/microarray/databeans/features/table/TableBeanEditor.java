package fi.csc.microarray.databeans.features.table;

import java.io.IOException;
import java.io.OutputStream;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.Session;
import fi.csc.microarray.databeans.Dataset;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.exception.MicroarrayException;

public class TableBeanEditor {
	
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(TableBeanEditor.class);

	
	private static final String COLUMNS_FEATURE = "/column/*";
	private Dataset bean;
	private EditableTable editableTable;

	public TableBeanEditor(Dataset bean) throws MicroarrayException {
		if (!bean.queryFeatures(COLUMNS_FEATURE).exists()) {
			throw new RuntimeException("bean " + bean.getName() + " is not tabular");
		}
		this.bean = bean;
		this.editableTable = toEditable();
		
	}
	
	public EditableTable getEditable() {
		return editableTable;
	}
	
	private EditableTable toEditable() throws MicroarrayException {
		Table table = bean.queryFeatures(COLUMNS_FEATURE).asTable();
		LinkedHashMap<String, List<String>> columns = new LinkedHashMap<String, List<String>>();  
		for (String columnName : table.getColumnNames()) {
			columns.put(columnName, new LinkedList<String>());
		}
		while (table.nextRow()) {
			for (String columnName : table.getColumnNames()) {
				String value = table.getStringValue(columnName);
				logger.debug("set editable column " + columnName + " to " + value);
				columns.get(columnName).add(value);
			}
		}
		
		// create editable table
		EditableTable editablePhenodataTable = new EditableTable();
		editablePhenodataTable.addColumns(columns);
		
		return editablePhenodataTable;
	}
	
	public void write() throws MicroarrayException, IOException {
		DataManager dataManager = Session.getSession().getDataManager();
		OutputStream out = dataManager.getContentOutputStreamAndLockDataBean(bean);
		try {
			getEditable().writeTo(out);
		} finally {
			dataManager.closeContentOutputStreamAndUnlockDataBean(bean, out);
		}
	}
}
