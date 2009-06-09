package fi.csc.microarray.client.operation.parameter;


/**
 * A parameter which takes integer values.
 * 
 * @author Janne KÃ¤ki
 *
 */
public class IntegerParameter extends Parameter {

	private int minValue;
	private int maxValue;
	private int value;
	
	/**
	 * Creates a new IntegerParameter with the given initial values.
	 * 
	 * @param name The name of the parameter.
	 * @param minValue The minimum allowed value of the parameter.
	 * @param maxValue The maximum allowed value of the parameter.
	 * @param initValue The initial value of the parameter.
	 * @throws IllegalArgumentException If the three initial values aren't
	 * 		   logically coherent (minimum value is bigger than maximum value,
	 * 		   or the init value is not between these limits).
	 */
	public IntegerParameter(
			String name, String description, int minValue, int maxValue, int initValue)
					throws IllegalArgumentException {
		super(name, description);
		this.minValue = minValue;
		if (maxValue < minValue) {
			throw new IllegalArgumentException("Minimum value for integer parameter " +
					this.getName() + " cannot be bigger than the maximum value.");
		}
		this.maxValue = maxValue;
		setIntegerValue(initValue);  // may throw IllegalArgumentException
	}
	
	public int getMinValue() {
		return minValue;
	}
	
	public int getMaxValue() {
		return maxValue;
	}
	
	public int getIntegerValue() {
		return value;
	}
	
	@Override
	public Object getValue() {
		return value;  // autowrapping really works, it seems :)
	}

	public void setMinValue(int newMinValue) {
		if (newMinValue > this.maxValue) {
			throw new IllegalArgumentException("New minimum value for " +
					this.getName() + " cannot exceed current maximum value.");
		}
		this.minValue = newMinValue;
		if (this.value < this.minValue) {
			this.value = this.minValue;
		}
	}
	
	public void setMaxValue(int newMaxValue) {
		if (newMaxValue < this.minValue) {
			throw new IllegalArgumentException("New maximum value for " +
					this.getName() + " cannot fall below current minimum value.");
		}
		this.maxValue = newMaxValue;
		if (this.value > this.maxValue) {
			this.value = this.maxValue;
		}
	}
	
	/**
	 * A method for setting the value as a pure integer, without the wrapping
	 * object around it.
	 * 
	 * @param newValue The new value for this parameter.
	 * @throws IllegalArgumentException If the suggested value was not within
	 * 		   the defined minimum and maximum limits.
	 */
	public void setIntegerValue(int newValue) throws IllegalArgumentException {
		if (newValue < minValue || newValue > maxValue) {
			throw new IllegalArgumentException("New value for integer parameter " +
					this.getName() + " must be inside given limits.");
		}
		this.value = newValue;
	}
	
	@Override
	public void setValue(Object newValue) throws IllegalArgumentException {
		if (newValue instanceof Integer) {
			setIntegerValue((Integer) newValue);
		} else {
			throw new IllegalArgumentException(newValue + " is an illegal " +
					"value for integer parameter " + this.getName() + ".");
		}
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

	@Override
	public String toString() {
		return this.getName() + ": " + value;
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
