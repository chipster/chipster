package fi.csc.microarray.client.cli;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;

import javax.jms.JMSException;

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
import fi.csc.microarray.messaging.auth.SimpleAuthenticationRequestListener;
import fi.csc.microarray.util.Strings;

public class CliClient {
	
	private static final String CONFIG = "config";
	private static final String USERNAME = "username";
	private static final String PASSWORD = "password";	
	
	private static final String VERBOSE = "verbose";
	private static final String QUIET = "quiet";
	private static final String YAML = "yaml";
	private static final String WORKING_COPY = "working-copy";
	
	private static final String INTERACTIVE = "interactive";
	private static final String EXIT = "exit";
	
	private static final String LIST = "list";
	private static final String VIEW = "view";
	private static final String PRINT = "print";
	private static final String HISTORY = "history";
	private static final String RENAME = "rename";
	private static final String DELETE = "delete";
	private static final String IMPORT = "import";
	private static final String EXPORT = "export";
	
	private static final String TOOLS = "tools";
	private static final String SEARCH_TERM = "search-term";
	private static final String TOOL = "tool";
	private static final String TOOL_ID = "tool-id";
	private static final String RUN = "run";
	private static final String DATASET = "dataset";
	private static final String FILE = "file";
	private static final String OLD_NAME = "old-name";
	private static final String NEW_NAME = "new-name";
	private static final String PARAMETER = "parameter";
	private static final String SESSION = "session";
	private static final String SAVE_WORKFLOW = "save-workflow";
	private static final String RUN_WORKFLOW = "run-workflow";
	
	private static final String OPEN_SESSION = "open-session";
	private static final String SAVE_SESSION = "save-session";
	private static final String LIST_SESSIONS = "list-sessions";
	private static final String DELETE_SESSION = "delete-session";
	private static final String CLEAR_SESSION = "clear-session";
		
	private static final String DEFAULT_WORKING_COPY = "cli-working-copy.zip";
	
	// without headless mode OSX will show this process in the dock and grab the focus
    static {
        System.setProperty("java.awt.headless", "true");
     }
	
	public static void main(String[] args) {
		new CliClient(args).runCliClient();
	}

	private Namespace nameSpace;
	private CliClientApplication app;
	private String[] args;
	private ArgumentParser parser;
	
	public CliClient(String[] args) {
		this.args = args;
	}
	
	private void runCliClient() {
		
		int exitValue = 1;
		try {
			this.parse();
			exitValue = 0;
		} catch (UserErrorException e) {
			System.err.println(e.getMessage());
		} catch (Exception e) {
			e.printStackTrace();
		}

		if (app != null) {
			app.quit();
		}
		System.exit(exitValue);
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
        
        addStringOption(parser, "-c", CONFIG,  "chipster client configuration file");
        addStringOption(parser, "-u", USERNAME, "chipster username");
        addStringOption(parser, "-p", PASSWORD, "chipster password");
        addStringOption(parser, "-W", WORKING_COPY, "name of the working copy session, either zip or cloud session").setDefault(DEFAULT_WORKING_COPY);        
        
        addBooleanOption(parser, "-v", VERBOSE, "more verbose output");
        addBooleanOption(parser, "-q", QUIET, "uppress status messages and print only requested data");
        addBooleanOption(parser, "-y", YAML, "output in yaml format for programmatical access");                         
        
        Subparsers subparsers = parser.addSubparsers();
        subparsers.title("commands");
                
        addCommand(subparsers, LIST, "list datasets");        
        addCommand(subparsers, VIEW, "view dataset details", DATASET);
        addCommand(subparsers, PRINT, "output dataset contents", DATASET);
        addCommand(subparsers, HISTORY, "view history, set verbose to view also source codes", DATASET);
        addCommand(subparsers, RENAME, "rename dataset", OLD_NAME).addArgument(NEW_NAME);
        addCommand(subparsers, DELETE, "delete dataset", DATASET);
        addCommand(subparsers, IMPORT, "import file", FILE);
        addCommand(subparsers, EXPORT, "export dataset to file", DATASET);                               
        
        addCommand(subparsers, TOOLS, "list tools, search term is optional").addArgument(SEARCH_TERM).nargs("?");        
        addCommand(subparsers, TOOL, "show tool details", TOOL_ID);
        
        Subparser run = addCommand(subparsers, RUN, "run tool");
        run.addArgument(TOOL_ID).required(true);
        // options instead of positional arguments because these are lists
        run.addArgument("--" + DATASET).nargs("*").help("input dataset(s) for a tool");
        run.addArgument("--" + PARAMETER).nargs("*").help("set parameters for a tool, e.g. parameter=VALUE");
        // default help has tool-id in the end, which won't work
        run.usage("run [-h] tool-id [--dataset [DATASET [DATASET ...]]] [--parameter [PARAMETER [PARAMETER ...]]]");
        
        Subparser saveWorkflow = addCommand(subparsers, SAVE_WORKFLOW, "save workflow");
        saveWorkflow.addArgument(FILE).help("save workflow to this file").required(true);
        saveWorkflow.addArgument(DATASET).help("start saving from this dataset").required(true);
        
        Subparser runWorkflow = addCommand(subparsers, RUN_WORKFLOW, "run workflow");
        runWorkflow.addArgument(FILE).help("run workflow of this file").required(true);
        runWorkflow.addArgument(DATASET).help("start running from this dataset").required(true);
        
        addCommand(subparsers, OPEN_SESSION, "open zip session or cloud session", SESSION);
        addCommand(subparsers, SAVE_SESSION, "save zip session or cloud session", SESSION);
        addCommand(subparsers, CLEAR_SESSION, "delete all datasets of the working copy session");
        addCommand(subparsers, LIST_SESSIONS, "list cloud sessions");
        addCommand(subparsers, DELETE_SESSION, "delete cloud session", SESSION);
        
        addCommand(subparsers, INTERACTIVE, "enter interactive mode");
        addCommand(subparsers, EXIT, "quit interactive mode").aliases("quit");

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
			
		if (isCommand(INTERACTIVE)) {
			
			Scanner scanner = new Scanner(System.in);
			
			try {
				printlnStatus("Chipster command line client in interactive mode, type '-h' for help or 'exit' to quit");
				while(!Thread.currentThread().isInterrupted()) {
					// process command line args on the first round
					try {
						parseArgs();
						
						if (isCommand(EXIT)) {
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
		return nameSpace.getBoolean(option);
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
		
		boolean yaml = isBooleanOption(YAML);
				
		if (isCommand(LIST)) {
			listDatasets(yaml);
		}
		
		if (isCommand(LIST_SESSIONS)) {
			listSessions(yaml);
		}
		
		if (isCommand(DELETE_SESSION)) {
			deleteSession(nameSpace.getString(DATASET));
		}
		
		if (isCommand(VIEW)) {
			String dataset = nameSpace.getString(DATASET);
			viewDataset(dataset, yaml);
		}
		
		if (isCommand(PRINT)) {
			String dataset = nameSpace.getString(DATASET);
			printDataset(dataset);
		}
		
		if (isCommand(TOOLS)) {
			tools(nameSpace.getString(SEARCH_TERM), yaml);
		}
		
		if (isCommand(TOOL)) {
			String tool = nameSpace.getString(TOOL_ID);			
			tool(tool, yaml);
		}
		
		if (isCommand(EXPORT)) {
			String dataset = nameSpace.getString(DATASET);
			exportDataset(dataset);
		}
		
		if (isCommand(IMPORT)) {
			String filename = nameSpace.getString(FILE);
			importDataset(filename);
		}
		
		if (isCommand(RENAME)) {									
			renameDataset(nameSpace.getString(OLD_NAME), nameSpace.getString(NEW_NAME));
		}
		
		if (isCommand(RUN)) {
			
			String tool = nameSpace.getString(TOOL_ID);
			List<String> datasets = nameSpace.<String> getList(DATASET);
			List<String> parameters = nameSpace.<String> getList(PARAMETER);
			
			run(tool, datasets, parameters);
		}
		
		if (isCommand(CLEAR_SESSION)) {
			clearSession();
		}
		
		if (isCommand(SAVE_WORKFLOW)) {
			String data = nameSpace.getString(DATASET);
			String file = nameSpace.getString(FILE);
			saveWorkflow(data, file);
		}
		
		if (isCommand(RUN_WORKFLOW)) {
			String data = nameSpace.getString(DATASET);
			String file = nameSpace.getString(FILE);
			runWorkflow(data, file);
		}
		
		if (isCommand(DELETE)) {
			deleteDataset(nameSpace.getString(DATASET));
		}
			
		if (isCommand(HISTORY)) {
			String dataset = nameSpace.getString(DATASET);
			historyOfDataset(dataset, yaml);			
		}
		
		if (isCommand(OPEN_SESSION)) {
			openSession(nameSpace.getString(SESSION));
		}
						
		if (isCommand(SAVE_SESSION)) {
			saveSession(nameSpace.getString(SESSION));
		}
	}
	
	private void initClient() throws UserErrorException, IOException,
	IllegalConfigurationException, MicroarrayException {
		if (!isStringOption(CONFIG)) {
			throw new UserErrorException("config not set");
		}

		if (!isStringOption(USERNAME)) {
			throw new UserErrorException("username not set");
		}

		if (!isStringOption(PASSWORD)) {
			throw new UserErrorException("password not set");
		}

		DirectoryLayout.initialiseClientLayout(nameSpace.getString(CONFIG));
		
		SimpleAuthenticationRequestListener auth = new SimpleAuthenticationRequestListener(nameSpace.getString(USERNAME), nameSpace.getString(PASSWORD));		
		app = new CliClientApplication(auth, isBooleanOption(VERBOSE), isBooleanOption(QUIET));

		app.initialiseApplication(true);
	}
	
	private String openWorkingCopySession() throws UserErrorException, JMSException, Exception {
		
		String sessionName = nameSpace.getString(WORKING_COPY);
		
		if (isLocalSession(sessionName)) {
			File session = new File(sessionName);
			if (session.exists()) {
				// dataless session
				app.loadSessionAndWait(session, null, true, false, false);
			}
		} else {
			try {			
				String sessionId = getSessionId(sessionName); // throws UserErrorException if not found
				app.loadSessionAndWait(null, sessionId, true, false, false);
			} catch (UserErrorException e) {
				// working copy doesn't exist, will be created when saved
			}
		}		
		
		return sessionName;
	}

	private void openSession(String sessionName) throws UserErrorException, JMSException, Exception {
						
		if (isLocalSession(sessionName)) {
			File session = new File(sessionName);
			if (session.exists()) {
				app.loadSessionAndWait(session, null, false, false, false);
			} else {
				throw new UserErrorException("session not found: " + sessionName);
			}
		} else {
			String sessionId = getSessionId(sessionName); // throws UserErrorException if not found
			app.loadSessionAndWait(null, sessionId, true, false, false);
		}
	}

	private void checkCloudConfiguration() throws UserErrorException {
		if (!app.areCloudSessionsEnabled()) {
			throw new UserErrorException("cloud sessions are disabled on this server, use local .zip sessions instead");
		}
	}

	private void deleteSession(String sessionName) throws JMSException, Exception {
		String sessionId = getSessionId(sessionName);
		app.removeRemoteSession(sessionId);
	}

	private void saveWorkingCopySession(String workingCopy) throws Exception {
		if (isLocalSession(workingCopy)) {
			app.getDataManager().saveLightweightSession(new File(workingCopy));
		} else {
			saveSession(workingCopy);
		}
	}
	
	private void saveSession(String workingCopy) throws Exception {
		if (isLocalSession(workingCopy)) {
			app.saveSessionAndWait(false, new File(workingCopy), null);
		} else {
			checkCloudConfiguration();
			app.saveSessionAndWait(true, null, workingCopy);
		}
	}

	private boolean isLocalSession(String sessionName) {
		return sessionName.endsWith(".zip");
	}

	private void listSessions(boolean yaml) throws JMSException, Exception {
		List<DbSession> sessions = app.listRemoteSessions();
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
	
	private void clearSession() {		
		app.clearSessionWithoutConfirming();
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
		
		System.out.println(app.getHistoryText(bean, true, true, true, true, isBooleanOption(VERBOSE), true, true));
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
		if (!isBooleanOption(QUIET)) {
			System.out.print(status);
		}
	}
	
	private void printlnStatus(String status) {
		if (!isBooleanOption(QUIET)) {
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
		
		map.put("tool", oper.getID());
		map.put("name", oper.getFullName());
		map.put("description", oper.getDescription());
		map.put("help", oper.getHelpURL());
		
		ArrayList<HashMap<String, String>> parameters = getParameters(oper);	
		
		map.put("parameters", parameters);
		
		if (yaml) {
			dumpYaml(map);
		} else {
		
			System.out.println(map.get("tool") + "\t\t" + map.get("name"));
			System.out.println();					
			wrapAndPrint(map.get("description").toString());
			System.out.println(map.get("help"));				
			System.out.println("PARAMETERS");
			System.out.println();

			for (HashMap<String, String> paramMap : parameters) {

				System.out.println(paramMap.get("parameter") + "\t\t" + paramMap.get("name"));
				System.out.println("\t" + paramMap.get("type"));
				System.out.println();
				wrapAndPrint(paramMap.get("description"));
				System.out.println();
			}
		}
	}
	
	private void wrapAndPrint(String text) {
		String wrapped = WordUtils.wrap(text, 60);
		wrapped = "    " + wrapped.replace("\n", "\n    ");
		System.out.println(wrapped);
	}

	private ArrayList<HashMap<String, String>> getParameters(OperationDefinition oper) {
		
		ArrayList<HashMap<String, String>> parameters = new ArrayList<>();
		for (Parameter parameter : oper.getParameters()) {
			HashMap<String, String> parameterMap = new HashMap<String, String>();
			
			// we should create yaml data structures for different parameter types
			// just a string for now
			String typeText = parameterToString(parameter);
			
			parameterMap.put("parameter", parameter.getID());
			parameterMap.put("name", parameter.getDisplayName());
			parameterMap.put("type", typeText);
			parameterMap.put("description", parameter.getDescription());
						
			parameters.add(parameterMap);
		}
		return parameters;
	}

	private ArrayList<HashMap<String, String>> getParameters(OperationRecord oper) {
		
		ArrayList<HashMap<String, String>> parameters = new ArrayList<>();
		for (ParameterRecord parameter : oper.getParameters()) {
			HashMap<String, String> parameterMap = new HashMap<String, String>();			
			
			parameterMap.put("parameter", parameter.getNameID().getID());
			parameterMap.put("name", parameter.getNameID().getDisplayName());
			parameterMap.put("description", parameter.getNameID().getDescription());
			parameterMap.put("value", parameter.getValue());
						
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
					yamlTool.put("tool", toolId);
					yamlTool.put("name", toolName);
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
						System.out.println(yamlModuleName + "\t\t" + yamlCategoryName + "\t\t" + yamlTool.get("tool") + "\t\t" + yamlTool.get("name"));
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
		map.put("name", bean.getName());
		map.put("date", new DateTime(bean.getDate()).toString());
		map.put("size", bean.getSize()); // Long
		map.put("notes", bean.getNotes());
		
		OperationRecord oper = bean.getOperationRecord();
		map.put("tool", oper.getNameID().getID());
		
		ArrayList<String> inputs = new ArrayList<>();	
		for (InputRecord input : oper.getInputs()) {
			inputs.add(input.getValue().getName());
		}		
		map.put("inputs", inputs);
				
		ArrayList<HashMap<String, String>> params = getParameters(oper);
		map.put("parameters", params);
				
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

			System.out.println("Dataset            " + map.get("name"));
			System.out.println("Date               " + map.get("date"));
			System.out.println("Size               " + map.get("size") + " bytes");
			System.out.println("Notes              " + map.get("notes"));
			System.out.println("Produced by tool   " + map.get("tool"));
			
			System.out.print(  "Using inputs       ");			
			System.out.println(Strings.delimit(inputs, " "));
			
			System.out.print(  "Parameters         ");
			for (HashMap<String, String> param : params) {
				System.out.print(param.get("parameter") + "=" + param.get("value") + " ");
			}
			System.out.println(); // new line after loop
			
			System.out.println();			
			System.out.println("Input of following tools and the datasets produced");
			for (String outputTool : outputTools.keySet()) {
				System.out.print(outputTool + "\t\t");
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

	private String parameterToString(Parameter parameter) {
		if (parameter instanceof IntegerParameter) {
			IntegerParameter number = (IntegerParameter) parameter;
			String str = "integer, default " + number.getValueAsString();
			if (number.getMinValue() != null || number.getMaxValue() != null) {
				str += " (" + number.getMinValue() + "-" + number.getMaxValue() + ")";
			}
			return str;
			
		} else if (parameter instanceof DecimalParameter) {
			DecimalParameter number = (DecimalParameter) parameter;
			String str = "decimal, default " + number.getValueAsString();
			if (number.getMinValue() != null || number.getMaxValue() != null) {
				str += " (" + number.getMinValue() + "-" + number.getMaxValue() + ")";
			}
			return str;
			
		} else if (parameter instanceof PercentageParameter) {
			PercentageParameter number = (PercentageParameter) parameter;
			String str = "percentage, default " + number.getValueAsString();
			if (number.getMinValue() != null || number.getMaxValue() != null) {
				str += " (" + number.getMinValue() + "-" + number.getMaxValue() + ")";
			}
			return str;
			
		} else if (parameter instanceof MetaColnameParameter) {
			
			return "phenodata column name, default " + parameter.getValueAsString();
			
		} else if (parameter instanceof EnumParameter) {
			
			//also DataSelectionParameter
			
		    EnumParameter enumParam = (EnumParameter) parameter;
		    if (enumParam.getMaxCount() > 1) {
		        // List with multiple selections
		        return "multiple selection:\n" + enumToString(enumParam);
		    } else {
		        // List with single selection
		    	return "single selection:\n" + enumToString(enumParam);
		    }
		} else if (parameter instanceof StringParameter) {
			return "string, default " + parameter.getValueAsString();
			
		} else {		
			throw new IllegalArgumentException("The given Parameter object, " + parameter.getID() + ", was not of recognized type!");
		}
	}
	
	private String enumToString(EnumParameter enumParam) {
		
		String str = "";
		
		for (SelectionOption opt : ((SelectionOption[])enumParam.getOptions())) {
			str += "\t\t" + opt.getValue() + "\t\t" + opt.toString(); 
			if (opt.getValue().equals(enumParam.getValueAsString())) {
				str += "\t(default)";
			}
			str += "\n";
		}
		return str;
	}
	
	private DataBean getDataset(String name) throws UserErrorException {
		DataBean bean = app.getDataManager().getDataBean(name);
		if (bean == null) {
			throw new UserErrorException("dataset not found: " + name);			
		}
		return bean;
	}
	
	private String getSessionId(String sessionName) throws JMSException, Exception {
		List<DbSession> sessions = app.listRemoteSessions();
		
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
