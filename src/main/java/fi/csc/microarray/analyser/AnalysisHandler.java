package fi.csc.microarray.analyser;

import fi.csc.microarray.messaging.message.JobMessage;

public interface AnalysisHandler {

	public AnalysisJob createAnalysisJob(JobMessage message, AnalysisDescription description, ResultCallback resultHandler);

	public AnalysisDescription handle(String sourceResourceName) throws AnalysisException;
	
	public boolean supports(String sourceResourceName);
	
	public boolean isUptodate(AnalysisDescription description);
}
