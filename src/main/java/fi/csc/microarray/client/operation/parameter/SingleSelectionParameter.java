package fi.csc.microarray.client.operation.parameter;

import org.apache.log4j.Logger;

/**
 * A parameter that has a defined set of possible values, out of which only
 * one at a time can be selected.
 * 
 * @author Janne Käki
 *
 */
public class SingleSelectionParameter extends Parameter {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger
			.getLogger(SingleSelectionParameter.class);

	private SelectionOption[] options;
	private int selectedIndex;
	
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
		
		public static SelectionOption[] convertStrings(String[] strings) {
			SelectionOption[] options = new SelectionOption[strings.length];
			for (int i = 0; i < strings.length; i++) {
				options[i] = new SelectionOption(strings[i], strings[i]);
			}
			return options;
		}
	}
	/**
	 * Creates a new SingleSelectionParameter with the given initial values.
	 * 
	 * @param name The name of this parameter.
	 * @param options The array of all possible value objects of this parameter.
	 * @param selectedIndex The index of the initially selected value object.
	 * @throws IllegalArgumentException If the options array was null or empty,
	 * 		   or if the given initial index was out of the array's bounds.
	 */
	public SingleSelectionParameter(
			String name, String description, SelectionOption[] options, int selectedIndex)
					throws IllegalArgumentException {
		super(name, description);
		if (options == null) {
			throw new IllegalArgumentException("Options array for single " +
					"selection parameter " + name + " may not be null!");
		}
		if (options.length == 0) {
			throw new IllegalArgumentException("Options array for single " +
					"selection parameter " + name + " may not be empty!");
		}
		this.options = options;
		if (selectedIndex < 0 || selectedIndex >= options.length) {
			throw new IllegalArgumentException("Given default selection index " +
					"for parameter " + name + " was out of array bounds!");
		}
		this.selectedIndex = selectedIndex;
	}
	
	public SingleSelectionParameter(String name, String description) {
		super(name, description);
		this.options = null;
		this.selectedIndex = -1;
	}
	
	/**
	 * @return The array containing all the possible values of this parameter.
	 */
	public Object[] getOptions() {
		assert(options != null);
		return options;
	}
	
	/**
	 * @return The index of the currently selected value object.
	 */
	public int getSelectedIndex() {
		assert(selectedIndex != -1);
		return selectedIndex;
	}
	
	@Override
	public Object getValue() {
		assert(selectedIndex != -1);
		assert(options != null);
		assert(options[selectedIndex] != null);
		logger.debug("returning value " + options[selectedIndex].getValue() + " from index " + selectedIndex);
		return options[selectedIndex].getValue();
	}

	/**
	 * Resets the set of possible values of this parameter and selects one
	 * of them as default.
	 * 
	 * @param newOptions An array of the possible value objects.
	 * @param newSelectedIndex The index of the new selected value. 
	 * @throws IllegalArgumentException If the suggested selection index was
	 * 		   out of the bounds of the given value array.
	 */
	public void setOptions(SelectionOption[] newOptions, int newSelectedIndex)
			throws IllegalArgumentException {
		this.options = newOptions;
		setSelectedIndex(newSelectedIndex);		
	}

	/**
	 * Sets which one of this parameter's possible values is selected.
	 * 
	 * @param newSelectedIndex The index of the new selected value.
	 * @throws IllegalArgumentException If no such index is found in the
	 * 		   set of values.
	 */
	public void setSelectedIndex(int newSelectedIndex) 
			throws IllegalArgumentException {
		if (newSelectedIndex < 0 || newSelectedIndex >= this.options.length) {
			throw new IllegalArgumentException("given index " + newSelectedIndex + 
					" for parameter " + this.getName() + " is out of bounds");
		}
		this.selectedIndex = newSelectedIndex;
	}
	
	@Override
	public void setValue(Object newValue) {
		logger.debug("new value is " + newValue);
		for (int i = 0; i < options.length; i++) {
			if (options[i].getValue().equals(newValue.toString())) {
				setSelectedIndex(i);
				logger.debug("new index is " + i);
				return;
			}
		}
		throw new IllegalArgumentException("illegal value for parameter " + this.getName() + ": " + newValue.toString());
	}

	@Override
	public boolean checkValidityOf(Object valueObject) {
		return valueObject instanceof SelectionOption;
	}
	
	public String toString() {
		return this.getName() + ": " + options[selectedIndex].getValue();
	}

	@Override
	public String getValueAsJava() {
		return "\"" + options[selectedIndex].getValue() + "\"";
	}
	
	@Override
	public void parseValue(String stringValue) throws IllegalArgumentException {
		setValue(stringValue); // no parsing needed
	}

}
