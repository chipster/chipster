package fi.csc.microarray.comp;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;

import org.apache.log4j.Logger;

import fi.csc.chipster.toolbox.ToolboxTool;
import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.messaging.message.GenericJobMessage;

/**
 * Abstract base class for any JobFactory that connects to external interpreter to run commands.
 *  
 * @author Aleksi Kallio
 *
 */
public abstract class InterpreterJobFactory implements JobFactory {
	
	/**
	 * Logger for this class
	 */
	static final Logger logger = Logger.getLogger(InterpreterJobFactory.class);

	protected String interpreterCommand;
	protected String toolPath;
	protected String externalToolPath;
	protected ProcessPool processPool;
	protected boolean isDisabled = false;


	public InterpreterJobFactory(HashMap<String, String> parameters) throws IOException {
		
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
		int poolSizeMin = 5;
		int poolSizeMax = 20;
		int poolTimeout = 360;
		int processUseCountMax = 10;
		int processLifetimeMax = 36_000;

		try {
			// TODO fix R specificity
			Configuration configuration = DirectoryLayout.getInstance().getConfiguration();

			poolSizeMin = configuration.getInt("comp", "r-process-pool-size-min");
			poolSizeMax = configuration.getInt("comp", "r-process-pool-size-max");
			poolTimeout = configuration.getInt("comp", "r-process-pool-timeout");
			processUseCountMax = configuration.getInt("comp", "r-process-pool-process-use-count-max");
			processLifetimeMax = configuration.getInt("comp", "r-process-pool-process-lifetime-max");
		} catch (IllegalStateException e) {
			// DirectoryLayout isn't configured in RestCompServer
			// use above hard coded values for now, because this class is still in the old chipster project
			// and we can't use the new Config class here		
			logger.info("process pool config missing, using hard coded defaults (" + e.getMessage() + ")");
		}
		
		try {
			processPool = new ProcessPool(new File(parameters.get("workDir")), interpreterCommand, poolSizeMin, poolSizeMax, 
				poolTimeout, processUseCountMax, processLifetimeMax);
		} catch (Exception e) {
			this.isDisabled = true;
		}
	}

	@Override
	public abstract CompJob createCompJob(GenericJobMessage message, ToolboxTool tool, ResultCallback resultHandler) throws CompException;

	protected abstract String getStringDelimeter();
	protected abstract String getVariableNameSeparator();

	protected ToolDescription createToolDescription(ToolboxTool tool) throws CompException {

		File moduleDir = new File(tool.getModule());
		
		// create description
		ToolDescription ad;
		ad = new ToolDescriptionGenerator().generate(tool.getSadlDescription());

		// add interpreter specific stuff to ToolDescription
		ad.setCommand(interpreterCommand);
		ad.setImplementation(tool.getSource()); // include headers
		ad.setSourceCode(tool.getSource());

		// tool and script locations and other variables
		int threadsMax = 2;
		try {
			threadsMax = DirectoryLayout.getInstance().getConfiguration().getInt("comp", "job-threads-max");
		} catch (Exception e) {
			logger.warn("could not read job-threads-max from configuration", e);
		}
		int memoryMax = 8192;
		try {
			memoryMax = DirectoryLayout.getInstance().getConfiguration().getInt("comp", "job-memory-max");
		} catch (Exception e) {
			logger.warn("could not read job-memory-max from configuration", e);
		}

		
//		File commonScriptDir = new File(moduleDir.getParentFile(), "common" + toolPath);
		File modulesRootDir;;
		try {
			modulesRootDir = DirectoryLayout.getInstance().getModulesDir();
		} catch (IOException e) {
			throw new CompException(e);
		}
		
		// that toolPath is more like runtime path, e.g python or R
		File commonScriptDir = new File(modulesRootDir, "common" + toolPath);
		
		String vns = getVariableNameSeparator();
		String sd = getStringDelimeter();
		
		ad.setInitialiser(
				"chipster" + vns + "tools" + vns + "path = " + sd + externalToolPath + sd + "\n" +
				"chipster" + vns + "common" + vns + "path = " + sd + commonScriptDir.getAbsolutePath() + sd + "\n" + 
				"chipster" + vns + "module" + vns + "path = " + sd + new File(modulesRootDir, moduleDir.getName()) + sd + "\n" + 
				"chipster" + vns + "threads" + vns + "max = " + sd + threadsMax + sd + "\n" +
				"chipster" + vns + "memory" + vns + "max = " + sd + memoryMax + sd + "\n");

		return ad;		
	}

	
	public boolean isDisabled() {
		return this.isDisabled;
	}

}
