package fi.csc.microarray.analyser.java;

import fi.csc.microarray.analyser.AnalysisDescription;
import fi.csc.microarray.analyser.AnalysisException;
import fi.csc.microarray.analyser.AnalysisHandler;
import fi.csc.microarray.analyser.AnalysisJob;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.messaging.message.JobMessage;

public class JavaAnalysisHandler implements AnalysisHandler {

	@SuppressWarnings(value="unchecked")
	public AnalysisJob createAnalysisJob(JobMessage message, AnalysisDescription description, ResultCallback resultHandler) {
		try {
			Class<? extends Object> jobClass = (Class<? extends Object>)description.getImplementation();
			JavaAnalysisJobBase analysisJob = (JavaAnalysisJobBase)jobClass.newInstance();
			analysisJob.construct(message, description, resultHandler);
			return analysisJob;
			
		} catch (Exception e) {
			throw new RuntimeException("internal error: type " + description.getImplementation().toString() + " could not be instantiated");
		}
	}

	public AnalysisDescription handle(String sourceResourceName) throws AnalysisException {
		
		// get the job class
		Class<? extends Object> jobClass = null;
		try { 
			 jobClass = Class.forName(sourceResourceName);
		} catch (ClassNotFoundException e) {
			throw new AnalysisException("Could not load job class: " + sourceResourceName);
		}
		
		
		assert(JavaAnalysisJobBase.class.isAssignableFrom(jobClass));
		JavaAnalysisJobBase jobInstance; 
		try {
			jobInstance = (JavaAnalysisJobBase)jobClass.newInstance();
		} catch (Exception e) {
			// should not happen
			throw new RuntimeException(e);
		}
		
		AnalysisDescription description = new AnalysisDescription(this);
		description.setImplementation(jobClass);
		description.setVVSADL(jobInstance.getVVSADL());
		description.setSourceResourceName(jobClass.getName());
		
		return description;
	}


	public boolean supports(String sourceResourceName) {
		// get the job class
		Class<? extends Object> jobClass = null;
		try { 
			 jobClass = Class.forName(sourceResourceName);
		} catch (ClassNotFoundException e) {
			return false;
		}
		
		return JavaAnalysisJobBase.class.isAssignableFrom(jobClass);
	}

	public boolean isUptodate(AnalysisDescription description) {
		return true;
	}

}
