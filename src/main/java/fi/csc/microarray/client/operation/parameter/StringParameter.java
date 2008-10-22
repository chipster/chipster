package fi.csc.microarray.client.operation.parameter;

public class StringParameter extends Parameter {

	private String value = "";
	
	protected StringParameter(String name, String description, String initValue) {
		super(name, description);
		if (initValue != null) {
			setValue(initValue);
		}
	}

	@Override
	public Object getValue() {
		return value;
	}

	@Override
	public void setValue(Object newValue) throws IllegalArgumentException {
		assert(newValue instanceof String);
		this.value = (String)newValue;
	}

	@Override
	public boolean checkValidityOf(Object valueObject) {
		return valueObject instanceof String;
	}

	@Override
	public String toString() {
		return value;
	}

	@Override
	public String getValueAsJava() {
		return "\"" + value + "\"";
	}
	
	@Override
	public void parseValue(String stringValue) throws IllegalArgumentException {
		setValue(stringValue); // no parsing needed
	}
}
