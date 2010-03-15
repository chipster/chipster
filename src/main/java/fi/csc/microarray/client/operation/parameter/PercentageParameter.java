package fi.csc.microarray.client.operation.parameter;

/**
 * A parameter which takes percentage values. For now, percentages are handled
 * as integers between 0 and 100, perhaps a better precision is needed?
 * 
 * @author Janne KÃ¤ki
 *
 */
public class PercentageParameter extends Parameter {

	private final Integer minValue;
	private final Integer maxValue;
	private Integer value;
	
	/**
	 * Creates a new PercentageParameter with the given name and init value,
	 * using 0 and 100 as the minimum and maximum limits.
	 * 
	 * @param name The name of the parameter.
	 * @param initValue The initial value of the parameter.
	 * @throws IllegalArgumentException If the initial value was not a valid
	 * 		   percentage (between 0 and 100).
	 */
	public PercentageParameter(String name, int initValue)
			throws IllegalArgumentException {
		this(name, name, 0, 100, initValue);
	}
	
	/**
	 * Creates a new PercentageParameter with the given initial values.
	 * 
	 * @param name The name of the parameter.
	 * @param minValue The minimum allowed value of the parameter.
	 * @param maxValue The maximum allowed value of the parameter.
	 * @param initValue The initial value of the parameter.
	 * @throws IllegalArgumentException If the three initial values aren't
	 * 		   logically coherent (minimum value is bigger than maximum value,
	 * 		   or the init value is not between these limits, or any of these
	 * 		   values is not a valid percentage - between 0 and 100, that is).
	 */
	public PercentageParameter(
			String name, String description, Integer minValue, Integer maxValue, Integer initValue)
					throws IllegalArgumentException {
		super(name, description);
		
		if (minValue < 0) {
			throw new IllegalArgumentException(
					"Minimum value for percentage parameter " + this.getName() +
					" cannot be less than zero.");
		}
		if (maxValue > 100) {
			throw new IllegalArgumentException(
					"Maximum value for percentage parameter " + this.getName() +
					" cannot be over 100.");
		}
		this.minValue = minValue;
		if (maxValue < minValue) {
			throw new IllegalArgumentException(
					"Minimum value for percentage parameter " + this.getName() +
					" cannot be bigger than the maximum value.");
		}
		this.maxValue = maxValue;
		setIntegerValue(initValue);  // may throw IllegalArgumentException
	}
	
	public Integer getMinValue() {
		return minValue;
	}
	
	public Integer getMaxValue() {
		return maxValue;
	}

	public Integer getIntegerValue() {
		return value;
	}
	
	@Override
	public Object getValue() {
		return value;  // autowrapping really works, it seems :)
	}
	
	/**
	 * A method for setting the value as a pure integer, without the wrapping
	 * object around it.
	 * 
	 * @param newValue The new value for this parameter.
	 * @throws IllegalArgumentException If the suggested value was not within
	 * 		   the defined minimum and maximum limits.
	 */
	public void setIntegerValue(int newValue) {
		if (newValue < minValue || newValue > maxValue) {
			throw new IllegalArgumentException("Initial value for percentage " +
					"parameter " + this.getName() + " must be within given limits.");
		}
		this.value = newValue;
	}

	@Override
	public void setValue(Object newValue) {
		if (checkValidityOf(newValue) == true) {
			value = (Integer) newValue;
			return;
		}
		throw new IllegalArgumentException(newValue + " is an illegal " +
				"value for percentage parameter \"" + this.getName() + "\" (" + this.minValue + " ... " + this.maxValue + ")");
	}

	@Override
	public boolean checkValidityOf(Object valueObject) {
		if (!(valueObject instanceof Integer)) {
			return false;			
		}
		Integer intValue = (Integer) valueObject;
		if (intValue >= this.minValue && intValue <= this.maxValue) {
			return true;
		} else {
			return false;
		}
	}

	public String toString() {
		return this.getName() + ": " + value + "%";
	}

	@Override
	public String getValueAsJava() {
		return "" + value;
	}
	
	@Override
	public void parseValue(String stringValue) throws IllegalArgumentException {
		try {
			setValue(Integer.parseInt(stringValue));
			
		} catch (NumberFormatException e) {
			throw new IllegalArgumentException("cannot parse String value \"" + stringValue + "\"");
		}
	}

}
