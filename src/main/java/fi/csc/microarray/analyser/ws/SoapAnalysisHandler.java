package fi.csc.microarray.analyser.ws;

import fi.csc.microarray.analyser.AnalysisDescription;
import fi.csc.microarray.analyser.AnalysisJob;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.analyser.java.JavaAnalysisHandler;
import fi.csc.microarray.messaging.message.JobMessage;

public class SoapAnalysisHandler extends JavaAnalysisHandler {

	
	public AnalysisJob createAnalysisJob(JobMessage message, AnalysisDescription description, ResultCallback resultHandler) {
		AnalysisJob analysisJob = super.createAnalysisJob(message, description, resultHandler);
		if (!(analysisJob instanceof SoapAnalysisJob)) {
			throw new IllegalArgumentException("unsupported job referenced: " + analysisJob.getClass().getSimpleName() + " is not supported");
		}
		return analysisJob;
	}
}
