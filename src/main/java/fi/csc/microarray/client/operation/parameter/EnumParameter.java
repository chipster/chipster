package fi.csc.microarray.client.operation.parameter;

import java.util.Arrays;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

/**
 * A parameter that has a defined set of possible values, out of which only
 * one at a time can be selected.
 * 
 * @author Janne KÃ¤ki, naktinis
 *
 */
public class EnumParameter extends Parameter {
    /**
     * Logger for this class
     */
    private static final Logger logger = Logger
            .getLogger(EnumParameter.class);

    private SelectionOption[] options;
    private int minCount = 1;
    private int maxCount = 1;
    
    // A list of selected options
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
            return name;
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
    public EnumParameter(String name, String description, SelectionOption[] options,
                         List<SelectionOption> defaultOptions, int minCount, int maxCount)
                    throws IllegalArgumentException {
        super(name, description);
        if (options == null) {
            throw new IllegalArgumentException("Options array for a " +
                    "selection parameter " + name + " may not be null!");
        }
        if (options.length == 0) {
            throw new IllegalArgumentException("Options array for a " +
                    "selection parameter " + name + " may not be empty!");
        }
        this.options = options;
        
        // Check if defaults are ok by finding an intersection
        HashSet<SelectionOption> opts = new HashSet<SelectionOption>(Arrays.asList(options));
        opts.retainAll(defaultOptions);
        if (opts.size() != defaultOptions.size()) {
            throw new IllegalArgumentException("Some given default options " +
                                               "were incorrect!");            
        }
        
        // Mark default values as selected
        this.selectedOptions.addAll(defaultOptions);
        
        setMinCount(minCount);
        setMaxCount(maxCount);
    }
    
    public EnumParameter(String name, String description) {
        super(name, description);
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
                    this.getName() + " cannot exceed current maximum value.");
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
                    this.getName() + " cannot fall below current minimum value.");
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
        if (newValue instanceof List<?>) {
            // Multiple selections
            List<String> checkedList = (List<String>) newValue;
            for (SelectionOption option : options) {
                if (checkedList.contains(option.getValue())) {
                    selectedOptions.add(option);
                    logger.debug("adding value " + option.getValue());
                }
            }
            return;
        } else if (newValue instanceof SelectionOption) {
            // Single selection
            SelectionOption optionValue = (SelectionOption) newValue;
            for (SelectionOption option : options) {
                if (option.getValue().equals(optionValue.getValue())) {
                    selectedOptions.add(option);
                    logger.debug("new value is " + option.getValue());
                    return;
                }
            }
        }
        throw new IllegalArgumentException("illegal value for parameter " + this.getName() + ": " + newValue.toString());
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
     * @return comma separated String.
     */
    @Override
    public Object getValue() {
        assert(options != null);
        
        String selected = "";
        for (SelectionOption option : selectedOptions) {
            selected += "," + option.getValue(); 
        }
        
        if (!selected.equals("")) {
            // The first character is a comma
            selected = selected.substring(1);
        }
        
        logger.debug("returning value " + selected);
        return selected;
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
    
    public String toString() {
        return this.getName() + ": " + getValue();
    }

    @Override
    public String getValueAsJava() {
        return "\"" + getValue() + "\"";
    }
    
    @Override
    public void parseValue(String stringValue) throws IllegalArgumentException {
        setValue(stringValue); // no parsing needed
    }

}
