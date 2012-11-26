package fi.csc.microarray.databeans.features;

import fi.csc.microarray.exception.MicroarrayException;

public interface Feature {
	
	public String getName();
	public boolean exists();
	public Iterable<Float> asFloats() throws MicroarrayException;
	public Float asFloat() throws MicroarrayException;
	
	/**
	 * Returns Table presentation of the feature, or null, if not supported. 
	 * @return Table instance or null
	 */
	public Table asTable() throws MicroarrayException;
	public Iterable<String> asStrings() throws MicroarrayException;
	public String asString() throws MicroarrayException;	
}
