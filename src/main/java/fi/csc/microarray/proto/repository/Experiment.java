package fi.csc.microarray.proto.repository;

public interface Experiment {
	
	public String getName();

	public String getDescription();
	
	public String getUniqueIdentifier();
	
	/**
	 * Returns null if feature is not supported by the datasource.
	 */
	public String getExperimentType();

	/**
	 * Returns null if feature is not supported by the datasource.
	 */
	public String getExperimentDesign();
	
	/**
	 * Returns null if feature is not supported by the datasource.
	 */
	public String getAuthor();
	
	/**
	 * Returns null if feature is not supported by the datasource.
	 */
	public String getPublication();
	
	public Array getArray(String name);
	
	public Iterable<String> getArrayNames();

}
