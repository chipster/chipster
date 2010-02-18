package fi.csc.microarray.analyser.r;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.AnalysisDescription;
import fi.csc.microarray.analyser.AnalysisDescriptionGenerator;
import fi.csc.microarray.analyser.AnalysisException;
import fi.csc.microarray.analyser.AnalysisHandler;
import fi.csc.microarray.analyser.AnalysisJob;
import fi.csc.microarray.analyser.ProcessPool;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.analyser.VVSADLTool;
import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.description.VVSADLParser.ParseException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.module.chipster.ChipsterVVSADLParser;

public class RAnalysisHandler implements AnalysisHandler {
	private static final String FILETYPE = ".r";

	/**
	 * Logger for this class
	 */
	static final Logger logger = Logger
			.getLogger(RAnalysisHandler.class);

	private String rCommand;
	private String toolPath;
	private String externalToolPath;
	private String customScriptsDirName;
	private ProcessPool processPool;
	private boolean isDisabled = false;
	
	
	public RAnalysisHandler(HashMap<String, String> parameters) throws IOException {
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
		
		String command = parameters.get("command");
		if (command == null || command.equals("")) {
			throw new IllegalArgumentException("Illegal command string: " + command);
		}
		
		String commandParameters = parameters.get("commandParameters");
		
		if (commandParameters != null) {
			command += " " + commandParameters;
		}
		this.rCommand = command;
		this.toolPath = parameters.get("toolPath");
		this.customScriptsDirName = configuration.getString("comp", "custom-scripts-dir");
	
		this.externalToolPath = parameters.get("externalToolPath");
		if (externalToolPath == null) {
			throw new RuntimeException("externalToolPath must be set in runtimes.xml");
		}		
		
		// initialize process pool
		int poolSizeMin = configuration.getInt("comp", "r-process-pool-size-min");
		int poolSizeMax = configuration.getInt("comp", "r-process-pool-size-max");
		int poolTimeout = configuration.getInt("comp", "r-process-pool-timeout");
		int processUseCountMax = configuration.getInt("comp", "r-process-pool-process-use-count-max");
		int processLifetimeMax = configuration.getInt("comp", "r-process-pool-process-lifetime-max");

		try {
			processPool = new ProcessPool(new File(parameters.get("workDir")), rCommand, poolSizeMin, poolSizeMax, 
				poolTimeout, processUseCountMax, processLifetimeMax);
		} catch (Exception e) {
			this.isDisabled = true;
		}
	}
	
	public AnalysisJob createAnalysisJob(JobMessage message, AnalysisDescription description, ResultCallback resultHandler) {
		RAnalysisJob analysisJob = new RAnalysisJob();
		analysisJob.construct(message, description, resultHandler);
		analysisJob.setProcessPool(this.processPool);
		return analysisJob;
	}


	public AnalysisDescription handle(String sourceResourceName) throws AnalysisException {
		
		InputStream scriptSource;
		
		String scriptPath = toolPath + File.separator + sourceResourceName;
		logger.debug("creating descriptions from " + scriptPath);
		
		// check for custom script file
		File scriptFile = new File(customScriptsDirName + File.separator + scriptPath);
		if (scriptFile.exists()) {
			FileInputStream customScriptSource;
			try {
				customScriptSource = new FileInputStream(scriptFile);
			} catch (FileNotFoundException fnfe) {
				logger.error("Could not load custom script: " + scriptFile, fnfe);
				throw new AnalysisException("Could not load custom script: " + scriptFile);
			}
			scriptSource = customScriptSource;
			logger.info("using custom-script for " + scriptPath);
		} else {
			scriptSource = this.getClass().getResourceAsStream(scriptPath);
		}
		
		if (scriptSource == null) {
			throw new AnalysisException(scriptPath + " not found");
		}
		
		// read the VVSADL from the comment block in the beginning of file
		// and the actual source code
		VVSADLTool.ParsedRScript parsedScript;
		try {
			parsedScript = new VVSADLTool().parseRScript(scriptSource);
		} catch (MicroarrayException e) {				
			throw new AnalysisException(e);
		}
		
		// parse VVSADL and create AnalysisDescription		
		AnalysisDescription ad;
		try {
			ad = new AnalysisDescriptionGenerator().generate(new ChipsterVVSADLParser().parse(parsedScript.VVSADL), this);
		} catch (ParseException e) {
			throw new AnalysisException(e);
		}
		ad.setVVSADL(parsedScript.VVSADL);

		// add R specific stuff to AnalysisDescription
		ad.setCommand(rCommand);
		ad.setImplementation(parsedScript.rSource); // include headers
		ad.setSourceCode(parsedScript.rSource);
		ad.setSourceResourceName(sourceResourceName);
		ad.setSourceResourceFullPath(scriptPath);
		ad.setInitialiser("chipster.tools.path = '" + externalToolPath + "'\n");
		
		return ad;
	}

	
	public boolean supports(String sourceResourceName) {
		logger.debug("do we support " + sourceResourceName + ": " + sourceResourceName.toLowerCase().endsWith(FILETYPE));
		return sourceResourceName.toLowerCase().endsWith(FILETYPE);
	}


	public boolean isUptodate(AnalysisDescription description) {
		File scriptFile = new File(customScriptsDirName + description.getSourceResourceFullPath());
		
		// custom script exists and is than description creation
		if (scriptFile.exists()) {
			if (scriptFile.lastModified() > description.getCreationTime()) {
				return false;
			}
		} 
		
		// custom script has been deleted
		else if (description.isUpdatedSinceStartup()) {
			return false;
		}
		return true;
	}

	public boolean isDisabled() {
		return this.isDisabled;
	}
}
