package fi.csc.microarray.comp.java;

import java.io.IOException;
import java.util.HashMap;

import org.apache.log4j.Logger;

import fi.csc.chipster.toolbox.ToolboxTool;
import fi.csc.microarray.comp.CompException;
import fi.csc.microarray.comp.CompJob;
import fi.csc.microarray.comp.JobFactory;
import fi.csc.microarray.comp.ResultCallback;
import fi.csc.microarray.comp.ToolDescription;
import fi.csc.microarray.comp.ToolDescriptionGenerator;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLGenerator;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.module.chipster.ChipsterSADLParser;

public class JavaJobFactory implements JobFactory {

	/**
	 * Logger for this class
	 */
	static final Logger logger = Logger.getLogger(JavaJobFactory.class);
	
	private HashMap<String, String> parameters;

	
	public JavaJobFactory(HashMap<String, String> parameters) throws IOException, IllegalConfigurationException {
		this.parameters = parameters;
	}

	@SuppressWarnings(value="unchecked")
	public CompJob createCompJob(JobMessage message, ToolboxTool tool, ResultCallback resultHandler) throws CompException {
		ToolDescription description = createToolDescription(tool);
		
		try {
			Class<? extends Object> jobClass = (Class<? extends Object>)description.getImplementation();
			JavaCompJobBase analysisJob = (JavaCompJobBase)jobClass.newInstance();
			analysisJob.construct(message, description, resultHandler);
			return analysisJob;
			
		} catch (Exception e) {
			throw new RuntimeException("internal error: type " + description.getImplementation().toString() + " could not be instantiated");
		}
	}
	
	public HashMap<String, String> getParameters() {
		return parameters;
	}

	private ToolDescription createToolDescription(ToolboxTool tool) throws CompException {
		
		// get the job class
		Class<? extends Object> jobClass = null;
		try { 
			 jobClass = Class.forName(tool.getResourceName());
		} catch (ClassNotFoundException e) {
			logger.error("could not load job class: " + tool.getResourceName());
			throw new CompException("could not load job class: " + tool.getResourceName());
		}
		
		
		// parse SADL		
		SADLDescription sadlDescription;
		try {
			sadlDescription = new ChipsterSADLParser().parse(tool.getSADL());
		} catch (ParseException e) {
			throw new CompException(e);
		}
		
		// create analysis description
		ToolDescription ad;
		ad = new ToolDescriptionGenerator().generate(sadlDescription);
		
		// SADL back to string
		SADLGenerator.generate(sadlDescription);
		ad.setSADL(SADLGenerator.generate(sadlDescription));
		
		ad.setImplementation(jobClass);
		ad.setCommand("java");
		ad.setSourceCode("Source code for this tool is available within Chipster source code.");
		
		return ad;
	}


	@Override
	public boolean isDisabled() {
		return false;
	}

}

