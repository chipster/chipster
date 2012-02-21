package fi.csc.microarray.analyser.java;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.ToolDescription;
import fi.csc.microarray.analyser.ToolDescriptionGenerator;
import fi.csc.microarray.analyser.AnalysisException;
import fi.csc.microarray.analyser.AnalysisHandler;
import fi.csc.microarray.analyser.AnalysisJob;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLGenerator;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.module.chipster.ChipsterSADLParser;

public class JavaAnalysisHandler implements AnalysisHandler {

	/**
	 * Logger for this class
	 */
	static final Logger logger = Logger.getLogger(JavaAnalysisHandler.class);
	
	private HashMap<String, String> parameters;

	
	public JavaAnalysisHandler(HashMap<String, String> parameters) throws IOException, IllegalConfigurationException {
		this.parameters = parameters;
	}

	@SuppressWarnings(value="unchecked")
	public AnalysisJob createAnalysisJob(JobMessage message, ToolDescription description, ResultCallback resultHandler) {
		try {
			Class<? extends Object> jobClass = (Class<? extends Object>)description.getImplementation();
			JavaAnalysisJobBase analysisJob = (JavaAnalysisJobBase)jobClass.newInstance();
			analysisJob.construct(message, description, resultHandler);
			return analysisJob;
			
		} catch (Exception e) {
			throw new RuntimeException("internal error: type " + description.getImplementation().toString() + " could not be instantiated");
		}
	}
	
	public HashMap<String, String> getParameters() {
		return parameters;
	}

	public ToolDescription handle(File moduleDir, String sourceResourceName, Map<String, String> params) throws AnalysisException {
		
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
		
		// parse SADL		
		SADLDescription sadlDescription;
		try {
			sadlDescription = new ChipsterSADLParser().parse(jobInstance.getSADL(),
			        sourceResourceName);
		} catch (ParseException e) {
			throw new AnalysisException(e);
		}
		
		// create analysis description
		ToolDescription ad;
		ad = new ToolDescriptionGenerator().generate(sadlDescription, this);
		
		// SADL back to string
		SADLGenerator.generate(sadlDescription);
		ad.setSADL(SADLGenerator.generate(sadlDescription));
		
		ad.setImplementation(jobClass);
		ad.setCommand("java");
		ad.setToolFile(new File(jobClass.getCanonicalName()));
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

	public boolean isUptodate(ToolDescription description) {
		return true;
	}

	public boolean isDisabled() {
		return false;
	}

}

