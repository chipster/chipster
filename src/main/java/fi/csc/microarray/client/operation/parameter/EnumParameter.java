package fi.csc.microarray.client.operation.parameter;

import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

import fi.csc.microarray.util.Strings;

/**
 * A parameter that has a defined set of possible values
 * 
 * @author Janne KÃ¤ki, Rimvydas Naktinis, Aleksi Kallio
 *
 */
public class EnumParameter extends Parameter {
    /**
     * Logger for this class
     */
    private static final Logger logger = Logger
            .getLogger(EnumParameter.class);

    private int minCount = 1;
    private int maxCount = 1;

    /**
     * All selectable values of the parameter.
     */
    private SelectionOption[] options;

    /**
     * Currently selected values (subset of all selectable values)
     * 
     * @see #options
     */
    private List<SelectionOption> selectedOptions = new LinkedList<SelectionOption>();
    
    public static class SelectionOption {
        private String name;
        private String value;
        
        public SelectionOption(String name, String value) {
            this.name = name;
            this.value = value;
        }
        
        public String getValue() {
            return value;
        }
        
        public String toString() {
            if (name != null && !name.isEmpty()) {
            	return name;
            } else {
            	return value;
            }
        }
        
        public static SelectionOption[] convertStrings(String[] titles, String[] values) {
            SelectionOption[] options = new SelectionOption[titles.length];
            for (int i = 0; i < titles.length; i++) {
                options[i] = new SelectionOption(titles[i], values[i]);
            }
            return options;
        }
    }
    
    /**
     * Creates a new EnumParameter with the given initial values.
     * 
     * @param name The name of this parameter.
     * @param options The array of all possible value objects of this parameter.
     * @param selectedIndex The index of the initially selected value object.
     * @param minCount The minimum number of values that have to be selected.
     * @param maxCount The maximum number of values that can be selected,
     *        value larger than 0 indicates that it is represented by a multi-select
     *        component.
     * @throws IllegalArgumentException If the options array was null or empty,
     *            or if the given initial index was out of the array's bounds.
     */
    public EnumParameter(String id, String displayName, String description, SelectionOption[] options,
                         List<SelectionOption> defaultOptions, int minCount, int maxCount)
                    throws IllegalArgumentException {
        super(id, displayName, description);
        if (options == null) {
            throw new IllegalArgumentException("Options array for a " +
                    "selection parameter " + displayName + " may not be null!");
        }
        if (options.length == 0) {
            throw new IllegalArgumentException("Options array for a " +
                    "selection parameter " + displayName + " may not be empty!");
        }
        this.options = options;

        if (!defaultOptions.isEmpty()) {

        	// Check if defaults are ok by finding an intersection
        	HashSet<SelectionOption> opts = new HashSet<SelectionOption>(Arrays.asList(options));
        	opts.retainAll(defaultOptions);
        	if (opts.size() != defaultOptions.size()) {
        		throw new IllegalArgumentException("Some given default options " +
        		"were incorrect!");            
        	}

        	// Mark default values as selected
        	this.selectedOptions.addAll(defaultOptions);
        }
        
        setMinCount(minCount);
        setMaxCount(maxCount);
    }
    
    public EnumParameter clone() {
        EnumParameter param = (EnumParameter)super.clone();
        param.selectedOptions = new LinkedList<SelectionOption>(selectedOptions);
        return param;
    }
    
    public EnumParameter(String id, String displayName, String description) {
        super(id, displayName, description);
        this.options = null;
    }
    
    /**
     * @return The array containing all the possible values of this parameter.
     */
    public Object[] getOptions() {
        assert(options != null);
        return options;
    }
    
    /**
     * Get the minimum number of values that have to be chosen for this
     * parameter.
     */
    public int getMinCount() {
        return minCount;
    }
    
    /**
     * Get the maximum number of values that can be chosen for this
     * parameter.
     */
    public int getMaxCount() {
        return maxCount;
    }
    
    /**
     * Set the minimum number of values that have to be chosen for this
     * parameter.
     * 
     * @param newMinCount
     */
    public void setMinCount(int newMinCount) {
        if (newMinCount > this.maxCount) {
            throw new IllegalArgumentException("New minimum value for " +
                    this.getID() + " cannot exceed current maximum value.");
        }
        this.minCount = newMinCount;
    }
    
    /**
     * Set the maximum number of values that can be chosen for this
     * parameter.
     * 
     * @param newMinValue
     */
    public void setMaxCount(int newMaxCount) {
        if (newMaxCount < this.minCount) {
            throw new IllegalArgumentException("New maximum value for " +
                    this.getID() + " cannot fall below current minimum value.");
        }
        this.maxCount = newMaxCount;
    }
    
    /**
     * Set a value for this parameter. Value can be either
     * a list of strings representing multiple option values
     * or a SelectionOption representing a single choice.
     * 
     * @throws IllegalArgumentException if given element is not
     * found.
     */
    @SuppressWarnings("unchecked")
    @Override
    public void setValue(Object newValue) {
        
        // Remove old selections
        selectedOptions.clear();
        
        // Add new selections
        if (newValue instanceof List<?> || newValue instanceof String[]) {
            // Multiple selections
            List<String> checkedList = (List<String>) (newValue instanceof List<?> ?
                    newValue : Arrays.asList(newValue));
            for (SelectionOption option : options) {
                if (checkedList.contains(option.getValue())) {
                    selectedOptions.add(option);
                    logger.debug("adding value " + option.getValue());
                }
            }
            return;
        } else if (newValue instanceof SelectionOption || newValue instanceof String) {
            // Single selection (preprocessed SelectionOption or raw String)
            String optionValue = newValue instanceof SelectionOption ?
                    ((SelectionOption) newValue).getValue() : (String) newValue;
            for (SelectionOption option : options) {
                if (option.getValue().equals(optionValue)) {
                    selectedOptions.add(option);
                    logger.debug("new value is " + option.getValue());
                    return;
                }
            }
        }
        throw new IllegalArgumentException("illegal value for parameter " + this.getID() + ": " + newValue.toString());
    }
    
    /**
     * Set selected options for this parameter.
     */
    public void setSelectedOptions(List<SelectionOption> newOptions) {
        selectedOptions = newOptions;
    }
    
    /**
     * @return currently selected options.
     */
    public List<SelectionOption> getSelectedOptions() {
        return selectedOptions;
    }

    /**
     * @return comma separated String or null if no values are selected
     */
    @Override
    public Object getValue() {
        assert(options != null);
        
        // gather selected option values
        LinkedList<String> selectedValues = new LinkedList<String>();
        for (SelectionOption selectedOption : selectedOptions) {
        	selectedValues.add(selectedOption.getValue());
        }
        
        // create comma separated string of selected values
        return Strings.delimit(selectedValues, ",");
    }   

    /**
     * Resets the set of possible values of this parameter and selects one
     * of them as default.
     * 
     * @param newOptions An array of the possible value objects.
     * @param newDefaults A list of selected values. 
     * @throws IllegalArgumentException If the suggested selection index was
     *            out of the bounds of the given value array.
     */
    public void setOptions(SelectionOption[] newOptions, List<SelectionOption> newDefaults)
            throws IllegalArgumentException {
        this.options = newOptions;
        setSelectedOptions(newDefaults);        
    }

    @Override
    public boolean checkValidityOf(Object valueObject) {
        return valueObject instanceof SelectionOption;
    }
    
    @Override
    public String toString() {
        return this.getID() + ": " + getValue();
    }

    @Override
    public String getValueAsJava() {
        return "\"" + getValue() + "\"";
    }
    
    /**
     * Takes a comma-separated value list and selects them
     * for this parameter.
     */
    @Override
    public void parseValue(String stringValue) throws IllegalArgumentException {
        // Try splitting stringValue and pass it to setValue
        String[] stringValues = stringValue.split(",");
        setValue(stringValues);
    }

	@Override
	public String getValueAsString() {
		return (String)getValue(); // getValue always returns String
	}

}
