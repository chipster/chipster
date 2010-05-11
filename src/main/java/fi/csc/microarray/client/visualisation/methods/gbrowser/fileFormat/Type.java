package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

public enum Type {
	
	STRING(String.class), 
	LONG(Long.class), 
	FLOAT(Float.class), 
	NEWLINE(String.class);

	private Class<?> javaType;

	private Type(Class<?> javaType) {
		this.javaType = javaType;
	}

	public Class<?> getJavaType() {
		return javaType;
	}
}