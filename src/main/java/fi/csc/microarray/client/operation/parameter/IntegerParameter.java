package fi.csc.microarray.client.operation.parameter;


/**
 * A parameter which takes integer values.
 * 
 * @author Janne KÃ¤ki
 *
 */
public class IntegerParameter extends Parameter {

	private Integer minValue;
	private Integer maxValue;
	private Integer value;
	
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
			String id, String displayName, String description, Integer minValue, Integer maxValue, Integer initValue)
					throws IllegalArgumentException {
		super(id, displayName, description);
		
        this.minValue = minValue;
        this.maxValue = maxValue;
		
        if (initValue == null) {
            this.value = null;
            return;
        }
	      
		if (maxValue < minValue) {
			throw new IllegalArgumentException("Minimum value for integer parameter " +
					this.getID() + " cannot be bigger than the maximum value.");
		}
		
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

	public void setMinValue(int newMinValue) {
		if (newMinValue > this.maxValue) {
			throw new IllegalArgumentException("New minimum value for " +
					this.getID() + " cannot exceed current maximum value.");
		}
		this.minValue = newMinValue;
		if (this.value < this.minValue) {
			this.value = this.minValue;
		}
	}
	
	public void setMaxValue(int newMaxValue) {
		if (newMaxValue < this.minValue) {
			throw new IllegalArgumentException("New maximum value for " +
					this.getID() + " cannot fall below current minimum value.");
		}
		this.maxValue = newMaxValue;
		if (this.value > this.maxValue) {
			this.value = this.maxValue;
		}
	}
	
	/**
	 * A method for setting the value as an Integer.
	 * 
	 * @param newValue The new value for this parameter.
	 * @throws IllegalArgumentException If the suggested value was not within
	 * 		   the defined minimum and maximum limits.
	 */
	public void setIntegerValue(Integer newValue) throws IllegalArgumentException {
		if (newValue != null && (newValue < minValue || newValue > maxValue)) {
			throw new IllegalArgumentException("New value for integer parameter " +
					this.getID() + " must be inside given limits.");
		}
		this.value = newValue;
	}
	
	@Override
	public void setValue(Object newValue) throws IllegalArgumentException {
		if (newValue instanceof Integer || newValue == null) {
		    setIntegerValue((Integer) newValue);
		} else {
			throw new IllegalArgumentException(newValue + " is an illegal " +
					"value for integer parameter " + this.getID() + ".");
		}
	}

	@Override
	public boolean checkValidityOf(Object valueObject) {
	    
	    // Allow null values for unfilled fields
	    if (valueObject == null) {
	        return true;
	    }
	    
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
		return this.getID() + ": " + value;
	}

	@Override
	public String getValueAsJava() {
		return "" + value;
	}
	
	@Override
	public void parseValue(String stringValue) throws IllegalArgumentException {
	       
        // Empty string means that no value is set
	    // This is possible for non-required parameters
        if (stringValue == null || stringValue.equals("")) {
            setValue(null);
            return;
        }
        
        // Otherwise string should represent an integer
		try {
			setValue(Integer.parseInt(stringValue));
		} catch (NumberFormatException e) {
			throw new IllegalArgumentException("cannot parse String value \"" + stringValue + "\"");
		}
	}

	@Override
	public String getValueAsString() {
	    return value != null ? value.toString() : "";
	}
}
