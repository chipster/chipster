package fi.csc.microarray.client.operation.parameter;

import java.util.LinkedList;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.exception.MicroarrayException;

/**
 * Parameter type for selecting some piece of data (column, dataset, metacolumn..).
 * Selection options are dependent on data binding.
 * 
 * @author Aleksi Kallio
 *
 */
public abstract class DataSelectionParameter extends SingleSelectionParameter {

	protected String initValue;
	
	public DataSelectionParameter(String name, String description, String initValue) {
		super(name, description);
		this.initValue = initValue;
	}
	
	/**
	 * Option "EMPTY" is always added as the last one.
	 */
	protected void loadOptionsFromColumnNames(DataBean data) throws MicroarrayException {
		int initIndex = 0;
		LinkedList<String> colNames = new LinkedList<String>();
		
		if(data != null){
			Table columns = data.queryFeatures("/column/*").asTable();
			for (String columnName : columns.getColumnNames()) {
				colNames.add(columnName);
			}
		}
		colNames.add("EMPTY");
		for (int i = 0; i < colNames.size(); i++) {
			if (colNames.get(i).equals(initValue)) {
				initIndex = i;
				break;
			}
		}
		String[] strings = (String[]) colNames.toArray(new String[colNames.size()]);
		SelectionOption[] optionObjects = SelectionOption.convertStrings(strings);
		setOptions(optionObjects, initIndex);
	}
}
