package fi.csc.microarray.analyser;

import java.io.File;
import java.util.Map;

import fi.csc.microarray.messaging.message.JobMessage;

public interface AnalysisHandler {

	/**
	 * Creates analysis job object using the analysis description and job message, i.e., instantiates the analysis
	 * job described by the description and parameterised by the parameters contained in job message. 
	 */
	public AnalysisJob createAnalysisJob(JobMessage message, AnalysisDescription description, ResultCallback resultHandler) throws AnalysisException;

	/**
	 * Creates analysis description from a String sourceResourceName. Interpretation of sourceResourceName
	 * is analysis handler specific. Analyser guarantees that this method is called only with supported
	 * sourceResourceName's.
	 * @param toolFilename 
	 * 
	 * @see #supports(String)
	 */
	public AnalysisDescription handle(File moduleDir, String toolFilename, Map<String, String> parameters) throws AnalysisException;
	
	/**
	 * Checks if given analysis description is up to date. Analyser guarantees that this method is called
	 * only with supported descriptions.
	 */
	public boolean isUptodate(AnalysisDescription description);

	/**
	 * Returns true if the handler is unable to create jobs. Handler is still able to 
	 * create descriptions.
	 * 
	 * @return
	 */
	public boolean isDisabled();
}
