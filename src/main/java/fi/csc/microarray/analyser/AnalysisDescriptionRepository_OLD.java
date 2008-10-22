package fi.csc.microarray.analyser;

public interface AnalysisDescriptionRepository_OLD {

	public abstract AnalysisDescription getDescription(String name);

	public abstract Iterable<AnalysisDescription> getDescriptions();

	/**
	 * Returns one huge VVSADL block that contains all loaded analysis 
	 * descriptions.
	 * @return huge block
	 */
	public abstract StringBuffer serialiseAsStringBuffer();

}