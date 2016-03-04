package fi.csc.microarray.client.cli;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.InvocationTargetException;
import java.net.MalformedURLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;

import javax.jms.JMSException;
import javax.swing.SwingUtilities;

import net.sourceforge.argparse4j.ArgumentParsers;
import net.sourceforge.argparse4j.impl.action.StoreTrueArgumentAction;
import net.sourceforge.argparse4j.inf.Argument;
import net.sourceforge.argparse4j.inf.ArgumentParser;
import net.sourceforge.argparse4j.inf.ArgumentParserException;
import net.sourceforge.argparse4j.inf.Namespace;
import net.sourceforge.argparse4j.inf.Subparser;
import net.sourceforge.argparse4j.inf.Subparsers;
import net.sourceforge.argparse4j.internal.HelpScreenException;

import org.apache.commons.io.IOUtils;
import org.apache.commons.lang3.StringUtils;
import org.apache.commons.lang3.text.WordUtils;
import org.joda.time.DateTime;
import org.yaml.snakeyaml.DumperOptions;
import org.yaml.snakeyaml.DumperOptions.FlowStyle;
import org.yaml.snakeyaml.Yaml;

import fi.csc.microarray.client.dataimport.ImportItem;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.operation.OperationRecord.InputRecord;
import fi.csc.microarray.client.operation.OperationRecord.ParameterRecord;
import fi.csc.microarray.client.operation.ToolCategory;
import fi.csc.microarray.client.operation.ToolModule;
import fi.csc.microarray.client.operation.parameter.DecimalParameter;
import fi.csc.microarray.client.operation.parameter.EnumParameter;
import fi.csc.microarray.client.operation.parameter.EnumParameter.SelectionOption;
import fi.csc.microarray.client.operation.parameter.IntegerParameter;
import fi.csc.microarray.client.operation.parameter.MetaColnameParameter;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.client.operation.parameter.PercentageParameter;
import fi.csc.microarray.client.operation.parameter.StringParameter;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.ContentType;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.DataNotAvailableHandling;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.ChecksumInputStream;
import fi.csc.microarray.filebroker.DbSession;
import fi.csc.microarray.filebroker.FileBrokerException;
import fi.csc.microarray.messaging.auth.SimpleAuthenticationRequestListener;
import fi.csc.microarray.util.Strings;

/**
 * Run in project root with command
 * 
 * java -cp bin/eclipse:ext/lib/* fi.csc.microarray.client.cli.CliClient
 * 
 * @author klemela
 *
 */
public class CliClient {
	
	private static final String KEY_INPUTS = "inputs";
	private static final String KEY_NOTES = "notes";
	private static final String KEY_SIZE = "size";
	private static final String KEY_VALUE = "value";
	private static final String KEY_DATE = "date";
	private static final String KEY_PARAMETERS = "parameters";
	private static final String KEY_TOOL = "tool";
	private static final String KEY_HELP = "help";
	private static final String KEY_DESCRIPTION = "description";
	private static final String KEY_NAME = "name";
	private static final String KEY_PARAM = "parameter";
	private static final String KEY_OPTIONS = "options";
	private static final String KEY_MAX = "max";
	private static final String KEY_MIN = "min";
	private static final String KEY_DEFAULT = "default";
	private static final String KEY_TYPE = "type";
	private static final String KEY_OPTION = "option";
	private static final String KEY_INTEGER = "INTEGER";
	private static final String KEY_DECIMAL = "DECIMAL";
	private static final String KEY_PERCENT = "PERCENT";
	private static final String KEY_STRING = "STRING";
	private static final String KEY_ENUM = "ENUM";
	private static final String KEY_METACOLUMN_SEL = "METACOLUMN_SEL ";
	
	
	private static final String OPT_CONFIG = "config";
	private static final String OPT_USERNAME = "username";
	private static final String OPT_PASSWORD = "password";	
	
	private static final String OPT_VERBOSE = "verbose";
	private static final String OPT_QUIET = "quiet";
	private static final String OPT_YAML = "yaml";
	private static final String OPT_WORKING_COPY = "working-copy";
	private static final String OPT_LOCAL = "local";
	private static final String OPT_CLOUD = "cloud";
	private static final String OPT_SHOW_ID = "show-id";
	
	private static final String CMD_INTERACTIVE = "interactive";
	private static final String CMD_EXIT = "exit";
	
	private static final String CMD_LIST_DATASETS = "list-datasets";
	private static final String CMD_DATASET = "dataset";
	private static final String CMD_PRINT = "print";
	private static final String CMD_HISTORY = "history";
	private static final String CMD_RENAME = "rename";
	private static final String CMD_DELETE = "delete";
	private static final String CMD_IMPORT = "import";
	private static final String CMD_EXPORT = "export";
	
	private static final String CMD_LIST_TOOLS = "list-tools";
	private static final String CMD_TOOL = "tool";
	private static final String CMD_RUN = "run";
	private static final String CMD_SAVE_WORKFLOW = "save-workflow";
	private static final String CMD_RUN_WORKFLOW = "run-workflow";
	
	private static final String ARG_SEARCH_TERM = "search-term";
	private static final String ARG_TOOL_ID = "tool-id";
	private static final String ARG_DATASET = "dataset";
	private static final String ARG_FILE = "file";
	private static final String ARG_OLD_NAME = "old-name";
	private static final String ARG_NEW_NAME = "new-name";
	private static final String ARG_PARAMETER = "parameter";
	private static final String ARG_SESSION = "session";
	
	private static final String CMD_OPEN_SESSION = "open-session";
	private static final String CMD_SAVE_SESSION = "save-session";
	private static final String CMD_LIST_SESSIONS = "list-sessions";
	private static final String CMD_DELETE_SESSION = "delete-session";
	private static final String CMD_CLEAR_SESSION = "clear-session";
		
	private static final String DEFAULT_WORKING_COPY = "cli-working-copy.zip";
	
	// without headless mode OSX will show this process in the dock and grab the focus
    static {
        System.setProperty("java.awt.headless", "true");
     }
	
	public static void main(String[] args) throws InvocationTargetException, InterruptedException {
		new CliClient(args).runCliClient();
	}

	private Namespace nameSpace;
	private CliClientApplication app;
	private String[] args;
	private ArgumentParser parser;
	
	public CliClient(String[] args) {
		this.args = args;
	}
	
	private void runCliClient() throws InvocationTargetException, InterruptedException {
		
		int exitValue = 1;
		try {
			this.parse();
			exitValue = 0;
		} catch (UserErrorException e) {
			System.err.println(e.getMessage());
		} catch (Exception e) {
			e.printStackTrace();
		}

		// shutdown in EDT to prevent EDT callbacks from running code after 
		// messaging is closed
		SwingUtilities.invokeAndWait(new ShutdownRunnable(exitValue));
	}
	
	public class ShutdownRunnable implements Runnable {			
		private int exitValue;

		public ShutdownRunnable(int exitValue) {
			this.exitValue = exitValue;
		}

		@Override
		public void run() {				
			if (app != null) {
				app.quit();
			}
			System.exit(exitValue);
		}
	}
	
	private void parse() throws JMSException, Exception {
		/*
		 * Argument '@file' will read parameters from the file. From argparse4j manual:
		 * 
		 * "The each line of the file is treated as one argument. Please be aware 
		 * that trailing empty lines or line with only white spaces are also considered 
		 * as arguments, although it is not readily noticeable to the user. The empty 
		 * line is treated as empty string."
		 */
        parser = ArgumentParsers.newArgumentParser("Chipster command line client", true, "-", "@");
        
        addStringOption(parser, "-c", OPT_CONFIG,  "chipster client configuration file");
        addStringOption(parser, "-u", OPT_USERNAME, "chipster username");
        addStringOption(parser, "-p", OPT_PASSWORD, "chipster password");
        addStringOption(parser, "-W", OPT_WORKING_COPY, "name of the working copy session, either zip or cloud session").setDefault(DEFAULT_WORKING_COPY);        
        
        addBooleanOption(parser, "-v", OPT_VERBOSE, "more verbose output");
        addBooleanOption(parser, "-q", OPT_QUIET, "uppress status messages and print only requested data");
        addBooleanOption(parser, "-y", OPT_YAML, "output in yaml format for programmatical access");                         
        
        Subparsers subparsers = parser.addSubparsers();
        subparsers.title("commands");
                
        addCommand(subparsers, CMD_LIST_DATASETS, "list datasets");        
        addCommand(subparsers, CMD_DATASET, "view dataset details", ARG_DATASET);
        addCommand(subparsers, CMD_PRINT, "output dataset contents", ARG_DATASET);
        addCommand(subparsers, CMD_HISTORY, "view history, set verbose to view also source codes", ARG_DATASET);
        addCommand(subparsers, CMD_RENAME, "rename dataset", ARG_OLD_NAME).addArgument(ARG_NEW_NAME);
        addCommand(subparsers, CMD_DELETE, "delete dataset", ARG_DATASET);
        addCommand(subparsers, CMD_IMPORT, "import file", ARG_FILE);
        addCommand(subparsers, CMD_EXPORT, "export dataset to file", ARG_DATASET);                               
        
        addCommand(subparsers, CMD_LIST_TOOLS, "list tools, search term is optional").addArgument(ARG_SEARCH_TERM).nargs("?");        
        addCommand(subparsers, CMD_TOOL, "show tool details, set verbose to view parameter help texts.", ARG_TOOL_ID);
        
        Subparser run = addCommand(subparsers, CMD_RUN, "run tool");
        run.addArgument(ARG_TOOL_ID).required(true);
        // options instead of positional arguments because these are lists
        run.addArgument("--" + ARG_DATASET).nargs("*").help("input dataset(s) for a tool");
        run.addArgument("--" + ARG_PARAMETER).nargs("*").help("set parameters for a tool, e.g. parameter=VALUE");
        // default help has tool-id in the end, which won't work
        run.usage("run [-h] tool-id [--dataset [DATASET [DATASET ...]]] [--parameter [PARAMETER [PARAMETER ...]]]");
        
        Subparser saveWorkflow = addCommand(subparsers, CMD_SAVE_WORKFLOW, "save workflow");
        saveWorkflow.addArgument(ARG_FILE).help("save workflow to this file").required(true);
        saveWorkflow.addArgument(ARG_DATASET).help("start saving from this dataset").required(true);
        
        Subparser runWorkflow = addCommand(subparsers, CMD_RUN_WORKFLOW, "run workflow");
        runWorkflow.addArgument(ARG_FILE).help("run workflow of this file").required(true);
        runWorkflow.addArgument(ARG_DATASET).help("start running from this dataset").required(true);
        
        addCommand(subparsers, CMD_OPEN_SESSION, "open zip session or cloud session", ARG_SESSION, OPT_CLOUD, OPT_LOCAL);
        addCommand(subparsers, CMD_SAVE_SESSION, "save zip session or cloud session", ARG_SESSION, OPT_CLOUD, OPT_LOCAL, OPT_SHOW_ID);
        addCommand(subparsers, CMD_CLEAR_SESSION, "delete all datasets of the working copy session");
        addCommand(subparsers, CMD_LIST_SESSIONS, "list cloud sessions");
        addCommand(subparsers, CMD_DELETE_SESSION, "delete cloud session", ARG_SESSION, OPT_CLOUD, OPT_LOCAL);
        
        addCommand(subparsers, CMD_INTERACTIVE, "enter interactive mode");
        addCommand(subparsers, CMD_EXIT, "quit interactive mode").aliases("quit");

        parser.epilog("use 'COMMAND -h' to show command arguments");                    		
			            	
		// parse the command line arguments
        try {
        	parseArgs();
        } catch (HelpScreenException e) {
			// parser has printed help text, we are done
			return;
		}

		initClient();			
		String workingCopy = openWorkingCopySession();
			
		if (isCommand(CMD_INTERACTIVE)) {
			
			Scanner scanner = new Scanner(System.in);
			
			try {
				printlnStatus("Chipster command line client in interactive mode, type '-h' for help or 'exit' to quit");
				while(!Thread.currentThread().isInterrupted()) {
					// process command line args on the first round
					try {
						parseArgs();
						
						if (isCommand(CMD_EXIT)) {
							break;
						}
						
						execute();
						saveWorkingCopySession(workingCopy);
						
					} catch (UserErrorException e) {
						System.err.println(e.getMessage());
					} catch (HelpScreenException e) {
						// skip execute() when parser shows a help text, because
						// it doesn't update nameSpace
					}
					
					printStatus(">>>");
					String lineString = scanner.nextLine();				
					
					args = Strings.splitConsideringQuotes(lineString, ' ').toArray(new String[0]);
				}
			} finally {
				scanner.close();
			}
		
		} else {
			execute();
			saveWorkingCopySession(workingCopy);
		}			
	}
	
	private boolean isBooleanOption(String option) {
		Boolean value = nameSpace.getBoolean(option);
		if (value != null) {
			return value;
		}
		return false;
	}
	
	private boolean isStringOption(String option) {
		return nameSpace.getString(option) != null;
	}
	
	private boolean isCommand(String cmd) {
		return nameSpace.getAttrs().containsKey(cmd);
	}

	private Argument addBooleanOption(ArgumentParser parser, String shortOption, String longOption, String help) {
		return addStringOption(parser, shortOption, longOption, help).action(new StoreTrueArgumentAction());
	}

	private Argument addStringOption(ArgumentParser parser, String shortOption,
			String longOption, String help) {
				
		return parser.addArgument(shortOption, "--" + longOption).dest(longOption).help(help);
	}
	
	private Subparser addCommand(Subparsers subparsers, String command,
			String help, String argument, String... options) {
		
		Subparser subparser = addCommand(subparsers, command, help);
		subparser.addArgument(argument);
		for (String option : options) {
			subparser.addArgument("--" + option).dest(option).action(new StoreTrueArgumentAction());
		}
		return subparser;
	}

	private Subparser addCommand(Subparsers subparsers, String command,
			String help, String argument) {
		
		Subparser subparser = addCommand(subparsers, command, help);
		subparser.addArgument(argument);
		return subparser;
	}

	private Subparser addCommand(Subparsers subparsers, String command,
			String help) {
		
		return subparsers.addParser(command).help(help).setDefault(command, true);
	}

	/**
	 * @throws HelpScreenException help text printed, nameSpace is not updated!
	 * @throws UserErrorException
	 */
	private void parseArgs() throws HelpScreenException, UserErrorException {
		
		try {
			nameSpace = parser.parseArgs(args);		
			
		} catch (ArgumentParserException e) {
			if (e instanceof HelpScreenException) {
				throw (HelpScreenException)e;
			} else {
				throw new UserErrorException(e.getMessage());
			}
		}
	}

	private void execute() throws JMSException, Exception {
		
		boolean yaml = isBooleanOption(OPT_YAML);
				
		if (isCommand(CMD_LIST_DATASETS)) {
			listDatasets(yaml);
			
		} else if (isCommand(CMD_LIST_SESSIONS)) {
			listSessions(yaml);
			
		} else if (isCommand(CMD_DELETE_SESSION)) {
			deleteSession(nameSpace.getString(ARG_DATASET));			
			
		} else if (isCommand(CMD_PRINT)) {
			String dataset = nameSpace.getString(ARG_DATASET);
			printDataset(dataset);
			
		} else if (isCommand(CMD_LIST_TOOLS)) {
			tools(nameSpace.getString(ARG_SEARCH_TERM), yaml);
			
		} else if (isCommand(CMD_TOOL)) {
			String tool = nameSpace.getString(ARG_TOOL_ID);			
			tool(tool, yaml);
			
		} else if (isCommand(CMD_EXPORT)) {
			String dataset = nameSpace.getString(ARG_DATASET);
			exportDataset(dataset);
			
		} else if (isCommand(CMD_IMPORT)) {
			String filename = nameSpace.getString(ARG_FILE);
			importDataset(filename);
			
		} else if (isCommand(CMD_RENAME)) {									
			renameDataset(nameSpace.getString(ARG_OLD_NAME), nameSpace.getString(ARG_NEW_NAME));
			
		} else 	if (isCommand(CMD_RUN)) {
			
			String tool = nameSpace.getString(ARG_TOOL_ID);
			List<String> datasets = nameSpace.<String> getList(ARG_DATASET);
			List<String> parameters = nameSpace.<String> getList(ARG_PARAMETER);			
			run(tool, datasets, parameters);
			
		} else if (isCommand(CMD_CLEAR_SESSION)) {
			clearSession();
			
		} else if (isCommand(CMD_SAVE_WORKFLOW)) {
			String data = nameSpace.getString(ARG_DATASET);
			String file = nameSpace.getString(ARG_FILE);
			saveWorkflow(data, file);
			
		} else if (isCommand(CMD_RUN_WORKFLOW)) {
			String data = nameSpace.getString(ARG_DATASET);
			String file = nameSpace.getString(ARG_FILE);
			runWorkflow(data, file);
			
		} else if (isCommand(CMD_DELETE)) {
			deleteDataset(nameSpace.getString(ARG_DATASET));
			
		} else if (isCommand(CMD_HISTORY)) {
			String dataset = nameSpace.getString(ARG_DATASET);
			historyOfDataset(dataset, yaml);
			
		} else if (isCommand(CMD_OPEN_SESSION)) {
			openSession(nameSpace.getString(ARG_SESSION));
			
		} else if (isCommand(CMD_SAVE_SESSION)) {
			saveSession(nameSpace.getString(ARG_SESSION));
			
		} else if (isCommand(CMD_DATASET)) {
			/* this must be after all command having "dataset" argument, because
			 * isCommand() can't tell the difference betweeen command and 
			 * argument.
			 */ 		
			String dataset = nameSpace.getString(ARG_DATASET);
			viewDataset(dataset, yaml);
		}
	}
	
	private void initClient() throws UserErrorException, IOException,
	IllegalConfigurationException, MicroarrayException {
		if (!isStringOption(OPT_CONFIG)) {
			throw new UserErrorException("config not set");
		}

		if (!isStringOption(OPT_USERNAME)) {
			throw new UserErrorException("username not set");
		}

		if (!isStringOption(OPT_PASSWORD)) {
			throw new UserErrorException("password not set");
		}

		DirectoryLayout.initialiseClientLayout(nameSpace.getString(OPT_CONFIG));
		
		SimpleAuthenticationRequestListener auth = new SimpleAuthenticationRequestListener(nameSpace.getString(OPT_USERNAME), nameSpace.getString(OPT_PASSWORD));		
		app = new CliClientApplication(auth, isBooleanOption(OPT_VERBOSE), isBooleanOption(OPT_QUIET));

		app.initialiseApplication(true);
	}	

	private String openWorkingCopySession() throws UserErrorException, JMSException, Exception {
		
		String sessionName = nameSpace.getString(OPT_WORKING_COPY);
		
		if (isLocalSession(sessionName, false)) {
			File session = new File(sessionName);
			if (session.exists()) {
				// dataless session
				app.getSessionManager().loadSessionAndWait(session, null, true, false, false);
			}
		} else {
			try {			
				String sessionId = getSessionId(sessionName); // throws UserErrorException if not found
				app.getSessionManager().loadSessionAndWait(null, sessionId, true, false, false);
			} catch (UserErrorException e) {
				// working copy doesn't exist, will be created when saved
			}
		}		
		
		return sessionName;
	}

	private void openSession(String sessionName) throws UserErrorException, JMSException, Exception {
						
		if (isLocalSession(sessionName, true)) {
			File session = new File(sessionName);
			if (session.exists()) {
				app.getSessionManager().loadSessionAndWait(session, null, false, false, false);
			} else {
				throw new UserErrorException("session not found: " + sessionName);
			}
		} else {
			String sessionId = getSessionId(sessionName); // throws UserErrorException if not found
			app.getSessionManager().loadSessionAndWait(null, sessionId, true, false, false);
		}
	}

	private void checkCloudConfiguration() throws UserErrorException {
		if (!app.getSessionManager().areCloudSessionsEnabled()) {
			throw new UserErrorException("cloud sessions are disabled on this server, use local .zip sessions instead");
		}
	}

	private void deleteSession(String sessionName) throws JMSException, Exception {
		if (isLocalSession(sessionName, true)) {
			File session = new File(sessionName);
			session.delete();
		} else {
			String sessionId = getSessionId(sessionName);
			app.getSessionManager().removeRemoteSession(sessionId);
		}
	}

	private void saveWorkingCopySession(String workingCopy) throws Exception {
		if (isLocalSession(workingCopy, false)) {
			app.getSessionManager().saveLightweightSession(new File(workingCopy));
		} else {
			saveSession(workingCopy);
		}
	}
	
	private void saveSession(String workingCopy) throws Exception {
		if (isLocalSession(workingCopy, true)) {
			app.getSessionManager().saveSessionAndWait(false, new File(workingCopy), null);
		} else {
			checkCloudConfiguration();
			app.getSessionManager().saveSessionAndWait(true, null, workingCopy);
			if (isBooleanOption(OPT_SHOW_ID)) {
				System.out.println(app.getSessionManager().getSessionId());
			}
		}
	}

	private boolean isLocalSession(String sessionName, boolean enableOptions) {
		if (enableOptions) {
			if (isBooleanOption(OPT_CLOUD)) {
				return false;
			} else if (isBooleanOption(OPT_LOCAL)) {
				return true;
			}
		}
		return sessionName.endsWith(".zip");
	}

	private void listSessions(boolean yaml) throws JMSException, Exception {
		List<DbSession> sessions = app.getSessionManager().listRemoteSessions();
		ArrayList<String> list = new ArrayList<String>();
		for (DbSession session : sessions) {
			list.add(session.getName());
		}
		print(list, yaml);
	}

	private void deleteDataset(String name) throws UserErrorException {
		DataBean bean = getDataset(name);
		app.deleteDatasWithoutConfirming(bean);
	}
	
	private void clearSession() throws MalformedURLException, FileBrokerException {		
		app.getSessionManager().clearSessionWithoutConfirming();
	}
	
	private void saveWorkflow(String dataset, String filename) throws IOException, UserErrorException {
		DataBean bean = getDataset(dataset);
		app.getSelectionManager().selectSingle(bean, this);
		app.saveWorkflow(new File(filename));
	}
	
	private void runWorkflow(String dataset, String filename) throws IOException, UserErrorException, InterruptedException {
		DataBean bean = getDataset(dataset);
		printlnStatus("Running workflow...");
		app.getSelectionManager().selectSingle(bean, this);				
		app.runWorkflowAndWait(new File(filename).toURI().toURL());
	}

	private void historyOfDataset(String dataset, boolean yaml) throws UserErrorException {
		DataBean bean = getDataset(dataset);
		
		if (yaml) {
			System.err.println("yaml output format isn't impelemented for history");
		}
		
		System.out.println(app.getHistoryText(bean, true, true, true, true, isBooleanOption(OPT_VERBOSE), true, true));
	}

	private void renameDataset(String oldName, String newName) throws UserErrorException {
		DataBean old = getDataset(oldName);
		app.renameDataItem(old, newName);
	}

	private void run(String toolId, List<String> datasets, List<String> parameters) throws MicroarrayException, UserErrorException {
				
		ArrayList<DataBean> inputs = new ArrayList<>();
		
		if (datasets != null) {
			for(String name : datasets) {
				inputs.add(getDataset(name));
			}
		}
		
		OperationDefinition tool = getTool(toolId);
		Operation operation = new Operation(tool, inputs.toArray(new DataBean[0]));
					
		if (parameters != null) {
			
			for (String param : parameters) {
				
				String[] nameAndValue = param.split("=");

				try {
					String name = nameAndValue[0];
					String value = nameAndValue[1];

					Object valueObj = getParameterValueObject(operation.getParameter(name), value);
					operation.setParameter(name, valueObj);
				} catch (Exception e) {
					throw new UserErrorException("illegal parameter: " + param + " (" + e.toString() + ")");
				}					
			}
		}
		
		Task task = app.executeOperation(operation);
		
		printStatus("Running...");
		
		try {
			while (app.getTaskExecutor().getTasks(true, true).contains(task)) {
				printStatus(".");
				Thread.sleep(1000);
			}
		} catch (InterruptedException e) {
		}
		printlnStatus("");
	}

	private void printStatus(String status) {
		if (!isBooleanOption(OPT_QUIET)) {
			System.out.print(status);
		}
	}
	
	private void printlnStatus(String status) {
		if (!isBooleanOption(OPT_QUIET)) {
			System.out.println(status);
		}
	}

	private Object getParameterValueObject(Parameter parameter, String value) {

		if (parameter instanceof IntegerParameter) {			
			//IntegerParameter integer = (IntegerParameter) parameter;
			return Integer.parseInt(value);

		} else if (parameter instanceof DecimalParameter) {
			//DecimalParameter number = (DecimalParameter) parameter;
			return Float.parseFloat(value);

		} else if (parameter instanceof PercentageParameter) {
			//PercentageParameter number = (PercentageParameter) parameter;
			return Integer.parseInt(value);

		} else if (parameter instanceof EnumParameter) {			
			// only single selection for now
			EnumParameter enumParam = (EnumParameter) parameter;			

			for (SelectionOption opt : ((SelectionOption[])enumParam.getOptions())) {				
				if (opt.getValue().equals(value)) {
					return opt;
				}
			}
		} else if (parameter instanceof StringParameter) {
			//StringParameter stringParam = (StringParameter) parameter;
			return value;

		} else {		
			throw new IllegalArgumentException("The given Parameter object, " + parameter.getID() + ", was not of recognized type!");
		}
		return null;
	}

	private void exportDataset(String dataset) throws UserErrorException {
		DataBean source = getDataset(dataset);
		File destination = new File(dataset);
		app.exportToFileAndWait(source, destination);		
	}
	
	private void importDataset(String filename) {
		File file = new File(filename);
		ContentType type = app.getDataManager().guessContentType(file);
		ImportItem item = new ImportItem(file, file.getName(), type);
		ArrayList<ImportItem> group = new ArrayList<>();
		group.add(item);
		app.importGroupAndWait(group, null);
	}

	private void tool(String tool, boolean yaml) throws UserErrorException {
		OperationDefinition oper = getTool(tool);
			
		HashMap<String, Object> map = new HashMap<String, Object>();
		
		map.put(KEY_TOOL, oper.getID());
		map.put(KEY_NAME, oper.getFullName());
		map.put(KEY_DESCRIPTION, oper.getDescription());
		map.put(KEY_HELP, oper.getHelpURL());
		
		ArrayList<HashMap<String, Object>> parameters = getParameters(oper);	
		
		map.put(KEY_PARAMETERS, parameters);
		
		if (yaml) {
			dumpYaml(map);
		} else {
		
			System.out.print(StringUtils.rightPad((String) map.get(KEY_TOOL), 50));
			System.out.print(StringUtils.rightPad((String) map.get(KEY_NAME), 50));
			System.out.println();
			System.out.println();					
			wrapAndPrint(map.get(KEY_DESCRIPTION).toString());
			System.out.println(map.get(KEY_HELP));				
			System.out.println("PARAMETERS");
			System.out.println();

			for (HashMap<String, Object> paramMap : parameters) {
				
				System.out.print(parameterToString(paramMap));
				if (isBooleanOption(OPT_VERBOSE)) {
					System.out.println();
					wrapAndPrint((String) paramMap.get(KEY_DESCRIPTION));
				}
				System.out.println();
			}
		}
	}
	
	private void wrapAndPrint(String text) {
		String wrapped = WordUtils.wrap(text, 60);
		wrapped = "    " + wrapped.replace("\n", "\n    ");
		System.out.println(wrapped);
	}

	private ArrayList<HashMap<String, Object>> getParameters(OperationDefinition oper) {
		
		ArrayList<HashMap<String, Object>> parameters = new ArrayList<>();
		for (Parameter parameter : oper.getParameters()) {
			HashMap<String, Object> parameterMap = new HashMap<>();
			
			
			parameterMap.put(KEY_PARAM, parameter.getID());
			parameterMap.put(KEY_NAME, parameter.getDisplayName());
			parameterMap.put(KEY_DESCRIPTION, parameter.getDescription());
			HashMap<String, Object> type = parameterToYaml(parameter);
			parameterMap.putAll(type);
						
			parameters.add(parameterMap);
		}
		return parameters;
	}

	private ArrayList<HashMap<String, String>> getParameters(OperationRecord oper) {
		
		ArrayList<HashMap<String, String>> parameters = new ArrayList<>();
		for (ParameterRecord parameter : oper.getParameters()) {
			HashMap<String, String> parameterMap = new HashMap<String, String>();			
			
			parameterMap.put(KEY_PARAM, parameter.getNameID().getID());
			parameterMap.put(KEY_NAME, parameter.getNameID().getDisplayName());
			parameterMap.put(KEY_DESCRIPTION, parameter.getNameID().getDescription());
			parameterMap.put(KEY_VALUE, parameter.getValue());
						
			parameters.add(parameterMap);
		}
		return parameters;
	}

	private void tools(String searchTerm, boolean yaml) {		
		
		HashMap<String, YamlModule> yamlModules = new HashMap<>();		
		for (ToolModule chipsterModule : app.getToolModules()) {
			
			YamlModule yamlModule = new YamlModule();			
			for (ToolCategory chipsterCategory : chipsterModule.getVisibleCategories()) {
				
				YamlCategory yamlCategory = new YamlCategory();				
				for (OperationDefinition chipsterTool : chipsterCategory.getToolList()) {
					YamlTool yamlTool = new YamlTool();
					
					String toolId = chipsterTool.getID();
					String toolName = chipsterTool.getDisplayName();
					yamlTool.put(KEY_TOOL, toolId);
					yamlTool.put(KEY_NAME, toolName);
					if (searchTerm == null || 
							toolId.toLowerCase().contains(searchTerm.toLowerCase()) || 
							toolName.toLowerCase().contains(searchTerm.toLowerCase())) {
						yamlCategory.add(yamlTool);
					}
				}
				yamlModule.put(chipsterCategory.getName(), yamlCategory);
				
			}
			yamlModules.put(chipsterModule.getModuleName(), yamlModule);
			
		}
		
		if (yaml) {
			dumpYaml(yamlModules);
		} else {
			for (String yamlModuleName : yamlModules.keySet()) {
				YamlModule yamlModule = yamlModules.get(yamlModuleName);
				for (String yamlCategoryName : yamlModule.keySet()) {
					YamlCategory yamlCategory = yamlModule.get(yamlCategoryName);
					for (YamlTool yamlTool : yamlCategory) {
						System.out.print(StringUtils.rightPad(yamlModuleName, 20));
						System.out.print(StringUtils.rightPad(yamlCategoryName, 40));
						System.out.print(StringUtils.rightPad(yamlTool.get(CMD_TOOL), 50));
						System.out.print(StringUtils.rightPad(yamlTool.get(KEY_NAME), 40));
						System.out.println();
					}
				}					
			}
		}
	}
	
	// unwind nested generics
	private static class YamlTool extends HashMap<String, String> {}
	private static class YamlCategory extends ArrayList<YamlTool> {}
	private static class YamlModule extends HashMap<String, YamlCategory> {}
	

	private void printDataset(String dataset)
			throws IOException, UserErrorException {
					
		DataBean bean = getDataset(dataset);			
		ChecksumInputStream stream = app.getDataManager().getContentStream(bean, DataNotAvailableHandling.EXCEPTION_ON_NA);
					
		try {
		    IOUtils.copy(stream, System.out);
		} finally {
			stream.close();
		}
	}

	private void viewDataset(String dataset, boolean yaml) throws UserErrorException {
					
		DataBean bean = getDataset(dataset);
		HashMap<String, Object> map = new HashMap<String, Object>();
		map.put(KEY_NAME, bean.getName());
		map.put(KEY_DATE, new DateTime(bean.getDate()).toString());
		map.put(KEY_SIZE, bean.getSize()); // Long
		map.put(KEY_NOTES, bean.getNotes());
		
		OperationRecord oper = bean.getOperationRecord();
		map.put(KEY_TOOL, oper.getNameID().getID());
		
		ArrayList<String> inputs = new ArrayList<>();	
		for (InputRecord input : oper.getInputRecords()) {
			
			String inputName = "";
			if (input.getValue() != null) {
				inputName = input.getValue().getName();
			} else {
				inputName = input.getDataId();
			}
			inputs.add(inputName);
		}		
		map.put(KEY_INPUTS, inputs);
				
		ArrayList<HashMap<String, String>> params = getParameters(oper);
		map.put(KEY_PARAMETERS, params);
				
		HashMap<String, ArrayList<String>> outputTools = new HashMap<>();
		
		// group output datasets by operation
		for (DataBean result : bean.getLinkSources(Link.derivationalTypes())) {
			String outputTool = result.getOperationRecord().getNameID().getID();
			
			if (!outputTools.containsKey(outputTool)) {
				outputTools.put(outputTool, new ArrayList<String>());
			}
									
			outputTools.get(outputTool).add(result.getName());
		}
		
		map.put("input of", outputTools);		
		
		if (yaml) {
			dumpYaml(map);
		} else {

			System.out.println("Dataset            " + map.get(KEY_NAME));
			System.out.println("Date               " + map.get(KEY_DATE));
			System.out.println("Size               " + map.get(KEY_SIZE) + " bytes");
			System.out.println("Notes              " + map.get(KEY_NOTES));
			System.out.println("Produced by tool   " + map.get(CMD_TOOL));
			
			System.out.print(  "Using inputs       ");			
			System.out.println(Strings.delimit(inputs, " "));
			
			System.out.print(  "Parameters         ");
			for (HashMap<String, String> param : params) {
				System.out.print(param.get(KEY_PARAM) + "=" + param.get(KEY_VALUE) + " ");
			}
			System.out.println(); // new line after loop
			
			System.out.println();			
			System.out.print(StringUtils.rightPad("INPUT OF", 50));
			System.out.print("OUTPUT DATASETS");
			System.out.println();
			
			for (String outputTool : outputTools.keySet()) {
				System.out.print(StringUtils.rightPad(outputTool, 50));
				ArrayList<String> outputDatasets = outputTools.get(outputTool);
				System.out.println(Strings.delimit(outputDatasets, " "));	
			}
		}
	}

	private void listDatasets(boolean yaml) {
		ArrayList<String> list = new ArrayList<String>();
		for (DataBean bean : app.getDataManager().databeans()) {
			list.add(bean.getName());
		}
		print(list, yaml);		
	}

	private void print(ArrayList<String> list, boolean yaml) {
		if (yaml) {
			dumpYaml(list);
		} else {
			for (String item : list) {
				System.out.println(item);
			}
		}
	}
	
	private void dumpYaml(Object yaml) {
		DumperOptions options = new DumperOptions();
		options.setDefaultFlowStyle(FlowStyle.BLOCK);
		Yaml yamlLib = new Yaml(options);
		System.out.print(yamlLib.dump(yaml));
	}

	@SuppressWarnings("unchecked")
	private String parameterToString(HashMap<String, Object> map) {
		
		String str = "";
		str += StringUtils.rightPad((String) map.get(KEY_PARAM), 30);
		
		String type = "";
		
		switch ((String)map.get(KEY_TYPE)) {
		case KEY_INTEGER:
			type = "integer" + minMaxDefaultToString(map);			
			break;
		case KEY_DECIMAL:
			type = "decimal" + minMaxDefaultToString(map);
			break;
		case KEY_PERCENT:
			type = "percent" + minMaxDefaultToString(map);
			break;
		case KEY_STRING:
			type = "string" + minMaxDefaultToString(map);
			break;
		case KEY_METACOLUMN_SEL:
			type = "phenodata column selection" + minMaxDefaultToString(map);			
			break;
		case KEY_ENUM:
			if (map.containsKey(KEY_MAX) && ((Integer)map.get(KEY_MAX)) > 1) {
				type = "multiple selection"; 				
			} else {
				type = "single selection";
			}
			if (map.containsKey(KEY_DEFAULT)) {
				type += ", default " + map.get(KEY_DEFAULT);
			}			
			break;
		}
		
		str += StringUtils.rightPad(type, 40);
		str += "\n";
		str += "    " + (String) map.get(KEY_NAME) + "\n";
		
		if (map.containsKey(KEY_OPTIONS)) {
			str += "\n";
			str += "    " + StringUtils.rightPad("OPTION", 60) + "NAME\n";
			for (HashMap<String, String> option: (ArrayList<HashMap<String, String>>)map.get(KEY_OPTIONS)) {
				str += "    ";
				str += StringUtils.rightPad(option.get(KEY_OPTION), 60);
				str += option.get(KEY_NAME) + "\n";
				//str += "  " + option.get(KEY_OPTION) + "\t\t\t\t\t\t\t" + option.get(KEY_NAME) + "\n";
			}
		}
		return str;
	}
	
	private String minMaxDefaultToString(HashMap<String, Object> map) {
		String str = "";
		if (map.containsKey(KEY_MIN) || map.containsKey(KEY_MAX)) {		
			str += " " + map.get(KEY_MIN) + "-" + map.get(KEY_MAX);
		}
		if (map.containsKey(KEY_DEFAULT)) {
			str += ", default " + map.get(KEY_DEFAULT);
		}
		return str;
	}

	private HashMap<String, Object> parameterToYaml(Parameter parameter) {
		HashMap<String, Object> map = new HashMap<>();
		if (parameter instanceof IntegerParameter) {
			map.put(KEY_TYPE, KEY_INTEGER);
			IntegerParameter number = (IntegerParameter) parameter;
			map.put(KEY_DEFAULT, number.getIntegerValue());
			if (number.getMinValue() != Integer.MIN_VALUE) {
				map.put(KEY_MIN, number.getMinValue());
			}
			if (number.getMaxValue() != Integer.MAX_VALUE) {
				map.put(KEY_MAX, number.getMaxValue());
			}
			return map;
			
		} else if (parameter instanceof DecimalParameter) {
			map.put(KEY_TYPE, KEY_DECIMAL);
			DecimalParameter number = (DecimalParameter) parameter;
			map.put(KEY_DEFAULT, number.getDecimalValue()); // Float
			if (number.getMinValue() != Float.MIN_VALUE) {
				map.put(KEY_MIN, number.getMinValue());
			}
			if (number.getMaxValue() != Float.MAX_VALUE) {
				map.put(KEY_MAX, number.getMaxValue());
			}
			return map;
			
		} else if (parameter instanceof PercentageParameter) {
			map.put(KEY_TYPE, KEY_PERCENT);
			PercentageParameter number = (PercentageParameter) parameter;
			map.put(KEY_DEFAULT, number.getIntegerValue());
			if (number.getMinValue() != Integer.MIN_VALUE) {
				map.put(KEY_MIN, number.getMinValue());
			}
			if (number.getMaxValue() != Integer.MAX_VALUE) {
				map.put(KEY_MAX, number.getMaxValue());
			}
			return map;
		} else if (parameter instanceof StringParameter) {
			map.put(KEY_TYPE, KEY_STRING);
			map.put(KEY_DEFAULT, parameter.getValueAsString());
			return map;
			
		// how about COLUMN_SEL?
			
		} else if (parameter instanceof MetaColnameParameter) {
			map.put(KEY_TYPE, KEY_METACOLUMN_SEL);
			map.put(KEY_DEFAULT, parameter.getValueAsString());
			return map;
			
		} else if (parameter instanceof EnumParameter) {
			
			//also DataSelectionParameter
			
		    EnumParameter enumParam = (EnumParameter) parameter;
		    map.put(KEY_TYPE, KEY_ENUM);		 
			map.put(KEY_DEFAULT, enumParam.getValueAsString());
			map.put(KEY_MIN, enumParam.getMinCount());
			map.put(KEY_MAX, enumParam.getMaxCount());
		    map.put(KEY_OPTIONS, enumOptionsToYaml(enumParam));
		    return map;
			
		} else {		
			throw new IllegalArgumentException("The given Parameter object, " + parameter.getID() + ", was not of recognized type!");
		}
	}
	
	private ArrayList<HashMap<String, String>> enumOptionsToYaml(EnumParameter enumParam) {
		ArrayList<HashMap<String, String>> options = new ArrayList<>();
		
		if (enumParam.getOptions() != null ) {
			for (SelectionOption opt : ((SelectionOption[])enumParam.getOptions())) {			
				HashMap<String, String> map = new HashMap<>();
				map.put(KEY_OPTION, opt.getValue());
				map.put(KEY_NAME, opt.toString());
				options.add(map);
			}
		}
		return options;
	}
	
	private DataBean getDataset(String name) throws UserErrorException {
		DataBean bean = app.getDataManager().getDataBean(name);
		if (bean == null) {
			throw new UserErrorException("dataset not found: " + name);			
		}
		return bean;
	}
	
	private String getSessionId(String sessionName) throws JMSException, Exception {
		List<DbSession> sessions = app.getSessionManager().listRemoteSessions();
		
		for (DbSession session : sessions) {
			if (sessionName.equals(session.getName())) {
				return session.getDataId();
			}
		}
		throw new UserErrorException("session not found: " + sessionName);		
	}
	
	private OperationDefinition getTool(String name) throws UserErrorException {
		OperationDefinition tool = app.getOperationDefinition(name);
		if (tool == null) {
			throw new UserErrorException("tool not found: " + name);
		}
		return tool;
	}
}
