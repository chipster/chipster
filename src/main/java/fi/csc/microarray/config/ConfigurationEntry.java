package fi.csc.microarray.config;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;


public class ConfigurationEntry {
	
	
	private static final String PROPERTY_REFERENCE_POSTFIX = "}";
	private static final String PROPERTY_REFERENCE_PREFIX = "${";

	enum Type {
		STRING,
		INT,
		BOOLEAN,
		STRINGS,
		INTS,
		BOOLEANS;

		public boolean isSingle() {
			return this == STRING || this == INT || this == BOOLEAN;
		}

		public static Type fromName(String typeName) {
			return valueOf(typeName.toUpperCase());
		}
	}
	
	private String[] stringValues;
	private int[] intValues;
	private boolean[] booleanValues;	
	private Type type;
	private String name;
	private boolean mustBeSet = false;
	
	public ConfigurationEntry(String name, String typeName) {
		this.type = Type.fromName(typeName);
		this.name = name;
	}
	
	public void setValue(String[] values) throws IllegalConfigurationException {
		
		// pre checks
		if (mustBeSet) {
			throw new IllegalConfigurationException(name + " not yet set");
		}
		
		if (type.isSingle() && values.length > 1) {
			throw new IllegalConfigurationException(name + " accepts only single value");
		}
		
		// do system property rewrites
		for (int i = 0; i < values.length; i++) {
			if (values[i].startsWith(PROPERTY_REFERENCE_PREFIX)) {
				String[] split = values[i].split(PROPERTY_REFERENCE_POSTFIX);
				String systemProperty = System.getProperty(split[0].substring(2));
				values[i] = systemProperty + split[1];
			}
		}
		
		// set value and check type
		switch (type) {
		case STRING:
		case STRINGS:
			this.stringValues = values;
			break;
			
		case INT:
		case INTS:
			this.intValues = new int[values.length];
			for (int i = 0; i < values.length; i++) {
				try {
					this.intValues[i] = Integer.parseInt(values[i]);
				} catch (NumberFormatException e) {
					throw new IllegalConfigurationException("illegal integer values " + values[i] + " for setting " + name);				
				}
			}
			break;
			
		case BOOLEAN:
		case BOOLEANS:
			this.booleanValues = new boolean[values.length];
			for (int i = 0; i < values.length; i++) {

				if (Boolean.TRUE.toString().equals(values[i].toLowerCase())) {
					this.booleanValues[i] = true;
				} else if (Boolean.FALSE.toString().equals(values[i].toLowerCase())) {
					this.booleanValues[i] = false;
				} else {
					throw new IllegalConfigurationException("illegal boolean value " + values[i] + " for setting " + name);
				}
			}

		default:
			throw new RuntimeException("unknown type " + type);	
		}
	}

	public String getString() {
		checkType(Type.STRING);

		if (stringValues.length > 1) {
			throw new IllegalArgumentException("multiple values found for " + name);

		} else {
			return stringValues[0];
		}
	}
	
	public int getInt() {
		checkType(Type.INT);
		
		if (intValues.length > 1) {
			throw new IllegalArgumentException("multiple values found for " + name);

		} else {
			return intValues[0];
		}
	}
	
	public boolean getBoolean() {
		checkType(Type.BOOLEAN);
		
		if (booleanValues.length > 1) {
			throw new IllegalArgumentException("multiple values found for " + name);

		} else {
			return booleanValues[0];
		}
	}
	
	public String[] getStrings() {
		checkType(Type.STRINGS);
		return stringValues;
	}
	
	public int[] getInts() {
		checkType(Type.INTS);
		return intValues;
	}
	
	public boolean[] getBooleans() {
		checkType(Type.BOOLEANS);
		return booleanValues;
	}
	
	private void checkType(Type suggested) {
		if (type != suggested) {
			throw new IllegalArgumentException(name + " has type " + type + ", not " + suggested);
			
		}
	}

	public void setMustBeSet(boolean mustBeSet) {
		this.mustBeSet = mustBeSet;		
	}
	
	public boolean mustBeSet() {
		return mustBeSet;
	}

	public String getName() {
		return name;
	}
}
