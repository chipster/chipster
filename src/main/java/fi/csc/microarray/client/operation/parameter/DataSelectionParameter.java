package fi.csc.microarray.client.operation.parameter;

import java.util.LinkedList;
import java.util.List;

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
public abstract class DataSelectionParameter extends EnumParameter {

	protected String initValue;
	
	public DataSelectionParameter(String id, String displayName, String description, String initValue) {
		super(id, displayName, description);
		this.initValue = initValue;
	}
	
	/**
	 * Option "EMPTY" is always added as the last one.
	 */
	protected void loadOptionsFromColumnNames(DataBean data) throws MicroarrayException {
		int initIndex = 0;
		LinkedList<String> colNames = new LinkedList<String>();
		
		if(data != null){
			try (Table columns = data.queryFeatures("/column/*").asTable()) {
				for (String columnName : columns.getColumnNames()) {
					colNames.add(columnName);
				}
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
		SelectionOption[] optionObjects = SelectionOption.convertStrings(strings, strings);
		
		// A list of default values
		List<SelectionOption> defaultOptions = new LinkedList<SelectionOption>();
		defaultOptions.add(optionObjects[initIndex]);
		
		setOptions(optionObjects, defaultOptions);
	}

    
    public void parseValueAndSetWithoutChecks(String stringValue) {
        // Try splitting stringValue and pass it to setValue
        String[] stringValues = stringValue.split(",");
        List<SelectionOption> selectedOptions = new LinkedList<SelectionOption>();
        for (String optionValue : stringValues) {
        	selectedOptions.add(new SelectionOption(optionValue, optionValue));
        }
        super.setSelectedOptions(selectedOptions);
    }


}
