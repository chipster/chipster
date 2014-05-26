package fi.csc.microarray.config;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;


public class ConfigurationEntry {
	
	private static final String PROPERTY_REFERENCE_POSTFIX = "}";
	private static final String PROPERTY_REFERENCE_PREFIX = "${";

	enum Type {
		STRING,
		INT,
		BOOLEAN;

		public static Type fromName(String typeName) {
			return valueOf(typeName.toUpperCase());
		}
	}
	
	private String stringValue;
	private int intValue;
	private boolean booleanValue;	
	private Type type;
	private String name;
	private boolean mustBeSet = false;
	
	public ConfigurationEntry(String name, String typeName) {
		this.type = Type.fromName(typeName);
		this.name = name;
	}
	
	public void setValue(String value) throws IllegalConfigurationException {
		
		// do system property rewrites
		if (value.startsWith(PROPERTY_REFERENCE_PREFIX)) {
			String[] split = value.split(PROPERTY_REFERENCE_POSTFIX);
			String systemProperty = System.getProperty(split[0].substring(2));
			value = systemProperty + split[1];
		}
		
		// set value and check type
		switch (type) {
		case STRING:
			this.stringValue = value;
			break;
			
		case INT:
			try {
				this.intValue = Integer.parseInt(value);
			} catch (NumberFormatException e) {
				throw new IllegalConfigurationException("illegal integer values " + value + " for setting " + name);				
			}
			break;
			
		case BOOLEAN:
			if (Boolean.TRUE.toString().equals(value.toLowerCase())) {
				this.booleanValue = true;
			} else if (Boolean.FALSE.toString().equals(value.toLowerCase())) {
				this.booleanValue = false;
			} else {
				throw new IllegalConfigurationException("illegal boolean value " + value + " for setting " + name);
			}
			break;

		default:
			throw new RuntimeException("unknown type " + type);	
		}
		
		// mark this value to be set
		this.mustBeSet = false;
	}

	public String getString() {
		checkType(Type.STRING);
		return stringValue;
	}
	
	public int getInt() {
		checkType(Type.INT);
		return intValue;
	}
	
	public boolean getBoolean() {
		checkType(Type.BOOLEAN);
		return booleanValue;
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
