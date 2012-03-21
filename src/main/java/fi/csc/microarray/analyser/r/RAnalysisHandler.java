package fi.csc.microarray.analyser.r;

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
import fi.csc.microarray.analyser.ProcessPool;
import fi.csc.microarray.analyser.ResultCallback;
import fi.csc.microarray.analyser.SADLTool;
import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.description.SADLDescription;
import fi.csc.microarray.description.SADLGenerator;
import fi.csc.microarray.description.SADLParser.ParseException;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.module.chipster.ChipsterSADLParser;

public class RAnalysisHandler implements AnalysisHandler {

	/**
	 * Logger for this class
	 */
	static final Logger logger = Logger
			.getLogger(RAnalysisHandler.class);

	private String rCommand;
	private String toolPath;
	private String externalToolPath;
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
	
	public AnalysisJob createAnalysisJob(JobMessage message, ToolDescription description, ResultCallback resultHandler) {
		RAnalysisJob analysisJob = new RAnalysisJob();
		analysisJob.construct(message, description, resultHandler);
		analysisJob.setProcessPool(this.processPool);
		return analysisJob;
	}


	public ToolDescription handle(File moduleDir, String toolFilename,
	                                  Map<String, String> params) throws AnalysisException {
		
		File toolFile = new File(moduleDir, toolPath + File.separator + toolFilename);
		
		InputStream scriptSource;
		
		try {
			scriptSource = new FileInputStream(toolFile);

		} catch (FileNotFoundException e) {
			toolFile = new File("modules/microarray", toolPath + File.separator + toolFilename);
			try {
				scriptSource = new FileInputStream(toolFile);

			} catch (FileNotFoundException fnfe) {
				throw new AnalysisException("script source " + toolFile + " not found.");
			}
		}
		
		// read the SADL from the comment block in the beginning of file
		// and the actual source code
		SADLTool.ParsedScript parsedScript;
		try {
			parsedScript = new SADLTool("#").parseScript(scriptSource);
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

		// add R specific stuff to ToolDescription
		ad.setCommand(rCommand);
		ad.setImplementation(parsedScript.source); // include headers
		ad.setSourceCode(parsedScript.source);
		ad.setToolFile(toolFile);
		
		// tool and script locations
		File commonRScriptDir = new File(moduleDir.getParentFile(), "common" + toolPath); 
		ad.setInitialiser("chipster.tools.path = '" + externalToolPath + "'\n" +
						  "chipster.common.path = '" + commonRScriptDir.getAbsolutePath() + "'\n" + 
						  "chipster.module.path = '" + moduleDir.getAbsolutePath() + "'\n");
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
