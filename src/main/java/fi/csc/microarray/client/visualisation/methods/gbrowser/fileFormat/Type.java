package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.text.DecimalFormat;

/**
 * Possible types of values stored in a file.
 *
 * @author Petri Klemelä
 */
public enum Type {
	
	STRING(String.class), 
	LONG(Long.class), 
	FLOAT(Float.class); 

	private static final DecimalFormat FLOAT_FORMAT = new DecimalFormat("0.#########");
	
	private Class<?> javaType;

	private Type(Class<?> javaType) {
		this.javaType = javaType;
	}

	public Class<?> getJavaType() {
		return javaType;
	}

	/**
	 * Does standardised conversion of different raw types to String
	 * Functionality has to be here because Java basic types cannot be
	 * extended for overriding toString.
	 */
	public static String toString(Object value) {
		if (value instanceof Float) {
			return FLOAT_FORMAT.format((Float)value);
			
		} else {
			return value.toString();
		}
	}
}