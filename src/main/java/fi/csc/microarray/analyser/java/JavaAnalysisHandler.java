package fi.csc.microarray.analyser.java;

import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.AnalysisDescription;
import fi.csc.microarray.analyser.AnalysisDescriptionGenerator;
import fi.csc.microarray.analyser.AnalysisException;
import fi.csc.microarray.analyser.AnalysisHandler;
import fi.csc.microarray.analyser.AnalysisJob;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.module.chipster.ChipsterSADLParser;

public class JavaAnalysisHandler implements AnalysisHandler {

	/**
	 * Logger for this class
	 */
	static final Logger logger = Logger.getLogger(JavaAnalysisHandler.class);

	
	public JavaAnalysisHandler(HashMap<String, String> parameters) throws IOException, IllegalConfigurationException {
		//Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
	}

	
	
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


	public AnalysisDescription handle(String sourceResourceName,
	                                  Map<String, String> params) throws AnalysisException {
		
		// get the job class
		Class<? extends Object> jobClass = null;
		try { 
			 jobClass = Class.forName(sourceResourceName);
		} catch (ClassNotFoundException e) {
			logger.error("could not load job class: " + sourceResourceName);
			throw new AnalysisException("could not load job class: " + sourceResourceName);
		}
		
		
		assert(JavaAnalysisJobBase.class.isAssignableFrom(jobClass));
		JavaAnalysisJobBase jobInstance; 
		try {
			jobInstance = (JavaAnalysisJobBase)jobClass.newInstance();
		} catch (Exception e) {
			logger.error("could not instantiate job: " + sourceResourceName);
			throw new RuntimeException(e);
		}
		
		// parse SADL and create AnalysisDescription		
		AnalysisDescription ad;
		try {
		    ad = new AnalysisDescriptionGenerator().generate(new ChipsterSADLParser().parse(jobInstance.getSADL()), this);
		} catch (ParseException e) {
			throw new AnalysisException(e);
		}
		
		ad.setImplementation(jobClass);
		ad.setCommand("java");
		ad.setSADL(jobInstance.getSADL());
		ad.setSourceResourceName(jobClass.getName());
		ad.setSourceResourceFullPath(jobClass.getCanonicalName());
		ad.setSourceCode("Source code for this tool is available within Chipster source code.");
		
		return ad;
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

	public boolean isDisabled() {
		return false;
	}

}

