package fi.csc.microarray.analyser.bsh;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.ToolDescription;
import fi.csc.microarray.analyser.ToolDescriptionGenerator;
import fi.csc.microarray.analyser.AnalysisException;
import fi.csc.microarray.analyser.AnalysisHandler;
import fi.csc.microarray.analyser.AnalysisJob;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.analyser.SADLTool;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLGenerator;
import fi.csc.microarray.description.SADLParser.ParseException;
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
	
	public BeanShellHandler(HashMap<String, String> parameters) throws IOException, IllegalConfigurationException {
		this.toolPath = parameters.get("toolPath");
	}
	
	public AnalysisJob createAnalysisJob(JobMessage message, ToolDescription description, ResultCallback resultHandler) {
		BeanShellJob analysisJob = new BeanShellJob();
		analysisJob.construct(message, description, resultHandler);
		return analysisJob;
	}


	public ToolDescription handle(File moduleDir, String toolFilename, Map<String, String> params) throws AnalysisException {
		
		File toolFile = new File(moduleDir, toolPath + File.separator + toolFilename);
		
		InputStream scriptSource;
		
		try {
			scriptSource = new FileInputStream(toolFile);
			
		} catch (FileNotFoundException e) {
			throw new AnalysisException("Script source " + toolFile + " not found.");
		}
		
		// read the SADL from the comment block in the beginning of file
		// and the actual source code
		SADLTool.ParsedScript parsedScript;
		try {
			parsedScript = new SADLTool("//").parseScript(scriptSource);
		} catch (IOException e) {				
			throw new AnalysisException(e);
		}
		// parse SADL		
		SADLDescription sadlDescription;
		try {
			sadlDescription = new ChipsterSADLParser().parse(parsedScript.SADL, toolFile.getName());
		} catch (ParseException e) {
			throw new AnalysisException(e);
		}
		
		// create analysis description
		ToolDescription ad;
		ad = new ToolDescriptionGenerator().generate(sadlDescription, this);
		
		// SADL back to string
		SADLGenerator.generate(sadlDescription);
		ad.setSADL(SADLGenerator.generate(sadlDescription));

		// add stuff to the ToolDescription
		ad.setCommand("BeanShell");
		ad.setImplementation(parsedScript.source); // include headers
		ad.setSourceCode(parsedScript.source);
//		ad.setSourceResourceName(sourceResourceName);
		ad.setSourceResourceFullPath(toolFile);
		
		return ad;
	}

	/**
	 * Check if the source file has been modified since the 
	 * ToolDescription was created.
	 */
	public boolean isUptodate(ToolDescription description) {
		File scriptFile = description.getToolFile();
		return scriptFile.lastModified() <= description.getCreationTime();
	}

	public boolean isDisabled() {
		return false;
	}
}
