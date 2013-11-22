package fi.csc.microarray.analyser;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLGenerator;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.module.chipster.ChipsterSADLParser;

/**
 * Abstract base class for any AnalysisHandler that connects to external interpreter to run commands.
 *  
 * @author Aleksi Kallio
 *
 */
public abstract class InterpreterAnalysisHandler implements AnalysisHandler {
	
	/**
	 * Logger for this class
	 */
	static final Logger logger = Logger.getLogger(InterpreterAnalysisHandler.class);

	protected String interpreterCommand;
	protected String toolPath;
	protected String externalToolPath;
	protected ProcessPool processPool;
	protected boolean isDisabled = false;


	@Override
	public abstract AnalysisJob createAnalysisJob(JobMessage message, ToolDescription description, ResultCallback resultHandler) throws AnalysisException;

	protected abstract String getSingleCommentLineSymbol();
	protected abstract String getStringDelimeter();
	protected abstract String getVariableNameSeparator();

	public InterpreterAnalysisHandler(HashMap<String, String> parameters) throws IOException {
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
		
		String command = parameters.get("command");
		if (command == null || command.equals("")) {
			throw new IllegalArgumentException("Illegal command string: " + command);
		}
		
		String commandParameters = parameters.get("commandParameters");
		
		if (commandParameters != null) {
			command += " " + commandParameters;
		}
		this.interpreterCommand = command;
		this.toolPath = parameters.get("toolPath");
	
		this.externalToolPath = parameters.get("externalToolPath");
		if (externalToolPath == null) {
			throw new RuntimeException("externalToolPath must be set in runtimes.xml");
		}		
		
		// initialize process pool
		// TODO fix R specificity
		int poolSizeMin = configuration.getInt("comp", "r-process-pool-size-min");
		int poolSizeMax = configuration.getInt("comp", "r-process-pool-size-max");
		int poolTimeout = configuration.getInt("comp", "r-process-pool-timeout");
		int processUseCountMax = configuration.getInt("comp", "r-process-pool-process-use-count-max");
		int processLifetimeMax = configuration.getInt("comp", "r-process-pool-process-lifetime-max");

		try {
			processPool = new ProcessPool(new File(parameters.get("workDir")), interpreterCommand, poolSizeMin, poolSizeMax, 
				poolTimeout, processUseCountMax, processLifetimeMax);
		} catch (Exception e) {
			this.isDisabled = true;
		}
	}
	

	@Override
	public ToolDescription handle(File moduleDir, String toolFilename,
			Map<String, String> params) throws AnalysisException {

		File toolFile = new File(moduleDir, toolPath + File.separator + toolFilename);

		InputStream scriptSource;

		try {
			scriptSource = new FileInputStream(toolFile);

		} catch (FileNotFoundException e) {
			throw new AnalysisException("script source " + toolFile + " not found.");
		}

		// read the SADL from the comment block in the beginning of file
		// and the actual source code
		SADLTool.ParsedScript parsedScript;
		try {
			parsedScript = new SADLTool(getSingleCommentLineSymbol()).parseScript(scriptSource);
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

		// add interpreter specific stuff to ToolDescription
		ad.setCommand(interpreterCommand);
		ad.setImplementation(parsedScript.source); // include headers
		ad.setSourceCode(parsedScript.source);
		ad.setToolFile(toolFile);
		ad.setModuleDir(moduleDir);

		// tool and script locations and other variables
		
		int threadsMax = 2;
		try {
			threadsMax = DirectoryLayout.getInstance().getConfiguration().getInt("comp", "job-threads-max");
		} catch (Exception e) {
			logger.warn("could not read job-threads-max from configuration", e);
		}
		File commonScriptDir = new File(moduleDir.getParentFile(), "common" + toolPath);
		
		String vns = getVariableNameSeparator();
		String sd = getStringDelimeter();
		
		ad.setInitialiser(
				"chipster" + vns + "tools" + vns + "path = " + sd + externalToolPath + sd + "\n" +
				"chipster" + vns + "common" + vns + "path = " + sd + commonScriptDir.getAbsolutePath() + sd + "\n" + 
				"chipster" + vns + "module" + vns + "path = " + sd + moduleDir.getAbsolutePath() + sd + "\n" + 
				"chipster" + vns + "threads" + vns + "max = " + sd + threadsMax + sd + "\n");
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
		return this.isDisabled;
	}

}
