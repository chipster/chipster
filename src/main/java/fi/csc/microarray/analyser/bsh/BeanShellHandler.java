package fi.csc.microarray.analyser.bsh;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.AnalysisDescription;
import fi.csc.microarray.analyser.AnalysisDescriptionGenerator;
import fi.csc.microarray.analyser.AnalysisException;
import fi.csc.microarray.analyser.AnalysisHandler;
import fi.csc.microarray.analyser.AnalysisJob;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.analyser.SADLTool;
import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.module.chipster.ChipsterSADLParser;
/**
 * Handler for the BeanShell jobs.
 * 
 * @author hupponen
 *
 */
public class BeanShellHandler implements AnalysisHandler {

	/**
	 * Logger for this class
	 */
	static final Logger logger = Logger.getLogger(BeanShellHandler.class);

	private final String toolPath;
	private final String customScriptsDirName;
	
	public BeanShellHandler(HashMap<String, String> parameters) throws IOException, IllegalConfigurationException {
		this.toolPath = parameters.get("toolPath");
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
		this.customScriptsDirName = configuration.getString("comp", "custom-scripts-dir");
	}
	
	public AnalysisJob createAnalysisJob(JobMessage message, AnalysisDescription description, ResultCallback resultHandler) {
		BeanShellJob analysisJob = new BeanShellJob();
		analysisJob.construct(message, description, resultHandler);
		return analysisJob;
	}


	public AnalysisDescription handle(String sourceResourceName,
	                                  Map<String, String> params) throws AnalysisException {
		
		InputStream scriptSource;
		
		String scriptPath = toolPath + File.separator + sourceResourceName;
		logger.debug("creating description from " + scriptPath);

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
			throw new AnalysisException("Script source " + sourceResourceName + " not found.");
		}
		
		// read the SADL from the comment block in the beginning of file
		// and the actual source code
		SADLTool.ParsedRScript parsedScript;
		try {
			parsedScript = new SADLTool().parseRScript(scriptSource, "//");
		} catch (MicroarrayException e) {				
			throw new AnalysisException(e);
		}
		
		// parse SADL and create AnalysisDescription		
		AnalysisDescription ad;
		try {
			ad = new AnalysisDescriptionGenerator().generate(new ChipsterSADLParser().parse(parsedScript.SADL), this);
		} catch (ParseException e) {
			throw new AnalysisException(e);
		}
		ad.setSADL(parsedScript.SADL);

		// add stuff to the AnalysisDescription
		ad.setCommand("BeanShell");
		ad.setImplementation(parsedScript.rSource); // include headers
		ad.setSourceCode(parsedScript.rSource);
		ad.setSourceResourceName(sourceResourceName);
		ad.setSourceResourceFullPath(scriptPath);
		
		return ad;
	}

	/**
	 * Check if the source file has been modified since the 
	 * AnalysisDescription was created.
	 */
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
		return false;
	}
}
