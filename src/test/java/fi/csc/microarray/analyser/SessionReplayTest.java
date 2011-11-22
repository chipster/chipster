package fi.csc.microarray.analyser;

import java.awt.event.MouseEvent;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URL;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.Map.Entry;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import javax.swing.Icon;

import org.springframework.validation.Errors;
import org.testng.Assert;

import fi.csc.microarray.analyser.AnalysisTestBase.JobResultListener;
import fi.csc.microarray.client.AtEndListener;
import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.RemoteServiceAccessor;
import fi.csc.microarray.client.ServiceAccessor;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dataimport.ImportItem;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.ChipsterDialog.PluginButton;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.operation.OperationRecord.InputRecord;
import fi.csc.microarray.client.operation.OperationRecord.ParameterRecord;
import fi.csc.microarray.client.operation.ToolCategory;
import fi.csc.microarray.client.operation.ToolModule;
import fi.csc.microarray.client.operation.parameter.DataSelectionParameter;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.tasks.Task.State;
import fi.csc.microarray.client.tasks.TaskException;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.MessagingTestBase;
import fi.csc.microarray.module.ModuleManager;
import fi.csc.microarray.module.chipster.MicroarrayModule;
import fi.csc.microarray.util.IOUtils;

public class SessionReplayTest extends MessagingTestBase {

	private static final String FLAG_FILE = "tool-test-ok";
	private static final String DEFAULT_SESSIONS_DIR = "sessions";
	
	private static final long TOOL_TEST_TIMEOUT = 1;
	private static final TimeUnit TOOL_TEST_TIMEOUT_UNIT = TimeUnit.HOURS;
	
	private static final boolean CHECK_EXACT_OUTPUT_SIZE = true;
	
	
	private List<ToolTestResult> toolTestResults = new LinkedList<ToolTestResult>();
	
	// Sessions which cause something to be thrown, normal tool failures etc not counted
	private LinkedHashMap<File, Throwable> sessionsWithErrors = new LinkedHashMap<File, Throwable>();
	
	private TaskExecutor executor;
	private DataManager manager, sourceManager;
	LinkedList<ToolModule> toolModules;
	
	public SessionReplayTest(String username, String password, String configURL) {
		super(username, password, configURL);
	}

	public boolean testSessions(String sessionsDirName) throws Exception {

		// Set up modules
		ModuleManager moduleManager = new ModuleManager("fi.csc.microarray.module.chipster.MicroarrayModule");
		Session.getSession().setModuleManager(moduleManager);
		
		// Set up main (target) system
		manager = new DataManager();
		moduleManager.plugAll(manager, null);
		executor = new TaskExecutor(super.endpoint, manager);
		toolModules = new LinkedList<ToolModule>();
		ServiceAccessor serviceAccessor = new RemoteServiceAccessor();
		serviceAccessor.initialise(manager, this.authenticationListener);
		serviceAccessor.fetchDescriptions(new MicroarrayModule());
		toolModules.addAll(serviceAccessor.getModules());
		Session.getSession().setClientApplication(new SessionLoadingSkeletonApplication(this, toolModules));
		
		// Set up source system
		sourceManager = new DataManager();
		moduleManager.plugAll(sourceManager, null);
		
		// Run all sessions
		File sessionsDir = new File(sessionsDirName);
		for (File testSession : sessionsDir.listFiles()) {
			
			// Zip files only, maybe should check if it really is a session file
			if (!testSession.getName().endsWith(".zip")) {
				continue;
			}

			// Clear data managers after previous session
			manager.deleteAllDataItems();
			sourceManager.deleteAllDataItems();
			
			// Test session
			try {
				testSession(testSession);
			} catch (Throwable e) {
				sessionsWithErrors.put(testSession, e);
			}
		}

		// Create reports
		createReports();
		
		// Get overall result
		boolean combinedResult = true;
		for (ToolTestResult toolTestResult : toolTestResults) {
			if (toolTestResult.getTestResult().equals(TestResult.FAIL)) {
				combinedResult = false;
			}
		}
		return combinedResult;
	}

	private void testSession(File session) throws IOException, MicroarrayException, TaskException, InterruptedException {

		Map<DataBean, DataBean> sourceDataBeanToTargetDataBean = new HashMap<DataBean, DataBean>();

		// Load session
		sourceManager.loadSession(session, false);
		
		// Pick import operations and copy imported data beans to target manager 
		// Also map OperationRecords to outputs TODO check that order is right, might need to traverse links
		LinkedList<OperationRecord> importOperationRecords = new LinkedList<OperationRecord>();
		Map<OperationRecord, List<DataBean>> outputMap = new HashMap<OperationRecord, List<DataBean>>();
		for (DataBean dataBean : sourceManager.databeans()) {
			OperationRecord operationRecord = dataBean.getOperationRecord();

			// pick import operations FIXME pick also any other without parent dataset
			if (OperationDefinition.IMPORT_DEFINITION_ID.equals(operationRecord.getNameID().getID())) {
				// copy imported databean, add mapping
				DataBean dataBeanCopy = manager.createDataBean(dataBean.getName(), session, dataBean.getContentUrl().getRef());
				sourceDataBeanToTargetDataBean.put(dataBean, dataBeanCopy);
				
				// avoid NPE 
				dataBeanCopy.setOperationRecord(operationRecord);
				
				// TODO what if not in the root folder in the source manager
				manager.getRootFolder().addChild(dataBeanCopy);
				importOperationRecords.add(operationRecord);
			}

			// store output mappings
			List<DataBean> outputs = outputMap.get(operationRecord);
			if (outputs != null) {
				outputs.add(dataBean);
			} else {
				outputs = new LinkedList<DataBean>();
				outputs.add(dataBean);
				outputMap.put(operationRecord, outputs);
			}
		}

		// Get operation records, avoid duplicates for tools with many outputs
		LinkedList<OperationRecord> operationRecords = new LinkedList<OperationRecord>();
		for (DataBean dataBean : sourceManager.databeans()) {
			OperationRecord operationRecord = dataBean.getOperationRecord();
			if (!operationRecords.contains(operationRecord)) {
				operationRecords.add(operationRecord);
			}
		}
	
		// Run operations in the order they were run originally
		for (OperationRecord operationRecord : operationRecords) {

			// Skip import operations
			if (importOperationRecords.contains(operationRecord)) {
				System.out.println("skipping import operation " + operationRecord.getFullName());
				continue;
			}

			// Get inputs
			LinkedList <DataBean> inputBeans = new LinkedList<DataBean>();
			for (InputRecord inputRecord : operationRecord.getInputs()) {
				DataBean inputBean = sourceDataBeanToTargetDataBean.get(inputRecord.getValue());
				if (inputBean != null) {
					inputBeans.add(inputBean);
				} else {
					String s = "not enough inputs for " + operationRecord.getFullName();
					System.out.println(s);
					throw new RuntimeException(s);
				}
			}

			// Set up task
			Operation operation = new Operation(getOperationDefinition(operationRecord.getNameID().getID(), toolModules), inputBeans.toArray(new DataBean[] {}));

			// Parameters, copy paste from workflows
			for (ParameterRecord parameterRecord : operationRecord.getParameters()) {
				if (parameterRecord.getValue() != null && !parameterRecord.getValue().equals("")) {	
					Parameter parameter = (Parameter)operation.getDefinition().getParameter(parameterRecord.getNameID().getID()).clone();
					if (parameter != null) {
						if (parameter instanceof DataSelectionParameter) {
							((DataSelectionParameter)parameter).parseValueAndSetWithoutChecks(parameterRecord.getValue());
						} else {
							parameter.parseValue(parameterRecord.getValue());
						}

						// Set it
						operation.setParameter(parameter.getID(), parameter.getValue());
					}
				}
			}

			Task task = executor.createTask(operation);
			
			// Execute the task
			System.out.println("running " + operation.getDefinition().getFullName());
			CountDownLatch latch = new CountDownLatch(1);
			task.addTaskEventListener(new JobResultListener(latch));
			executor.startExecuting(task);
			latch.await(TOOL_TEST_TIMEOUT, TOOL_TEST_TIMEOUT_UNIT);
			
			// Task failed
			if (!task.getState().equals(State.COMPLETED)) {
				
				// try to cancel if still running
				if (task.getState().equals(State.RUNNING)) {
					executor.kill(task);
					toolTestResults.add(new ToolTestResult(TestResult.FAIL, operation, session, task, "task did not finish before test timeout " +TOOL_TEST_TIMEOUT + " " + TOOL_TEST_TIMEOUT_UNIT.toString()));
					return;
				}
				toolTestResults.add(new ToolTestResult(TestResult.FAIL, operation, session, task, "task failed"));
				return;
			}
			
			// Check results
			// Target manager needs to be available through session for some of these to work
			Session.getSession().setDataManager(manager);
			try {
				
				// Link result beans, add to folders etc
				Session.getSession().getApplication().onFinishedTask(task, operation);

				// Check that number of results and result names match
				Iterator<DataBean> targetIterator = task.outputs().iterator();
				for (DataBean sourceBean : outputMap.get(operationRecord)) {
					if (targetIterator.hasNext()) {
						DataBean targetBean = targetIterator.next();
						if (!sourceBean.getName().equals(targetBean.getName())) {
							toolTestResults.add(new ToolTestResult(TestResult.FAIL, operation, session, task, "mismatch in result dataset names, "
									+ "expecting: " + sourceBean.getName() + " got: " + targetBean.getName()));
							return;
						}
					} 
					// Not enough results
					else {
						toolTestResults.add(new ToolTestResult(TestResult.FAIL, operation, session, task, "not enough result datasets"));
						return;
					}
				}
				if (targetIterator.hasNext()) {
					toolTestResults.add(new ToolTestResult(TestResult.FAIL, operation, session, task, "too many result datasets"));
					return;
				}

				// Find and replace metadata 
				targetIterator = task.outputs().iterator();
				for (DataBean sourceBean : outputMap.get(operationRecord)) {
					DataBean targetBean = targetIterator.next();

					// replace metadata contents from the source session
					if (Session.getSession().getPrimaryModule().isMetadata(targetBean)) {
						System.out.println("copying metadata for: " + targetBean.getName());
						OutputStream metadataOut = manager.getContentOutputStreamAndLockDataBean(targetBean);
						InputStream sourceIn = null;
						try {
							sourceIn = sourceBean.getContentByteStream();
							IOUtils.copy(sourceIn, metadataOut);
						} finally {
							IOUtils.closeIfPossible(sourceIn);
							manager.closeContentOutputStreamAndUnlockDataBean(targetBean, metadataOut);
						}
					}
				}

				// Compare data beans, add source bean -> target bean mapping
				targetIterator = task.outputs().iterator();
				for (DataBean sourceBean : outputMap.get(operationRecord)) {
					DataBean targetBean = targetIterator.next();
					try {
						compareDataBeans(sourceBean, targetBean);
					} catch (Throwable t) {
						toolTestResults.add(new ToolTestResult(TestResult.FAIL, operation, session, task, t.getMessage()));
						return;
					}
					
					// Add source bean -> target bean mapping, needed for further operations
					sourceDataBeanToTargetDataBean.put(sourceBean, targetBean);
					
					// Add result
					toolTestResults.add(new ToolTestResult(TestResult.OK, operation, session, task));
				
				}
			} finally {
				// Set session data manager back to null to avoid problems
				Session.getSession().setDataManager(null);
			}
		}
	}


	private OperationDefinition getOperationDefinition(String toolId, LinkedList<ToolModule> toolModules) {
		for (ToolModule module : toolModules) {
			OperationDefinition tool = module.getOperationDefinition(toolId);
			if (tool != null) {
				return tool;
			}
		}
		return null;
	}

	/**
	 * 
	 * @param bean1
	 * @param bean2
	 * 
	 * @throws Errors if not equal
	 */
	private void compareDataBeans(DataBean bean1, DataBean bean2) {

		// Check name, may be already checked though
		Assert.assertEquals(bean1.getName(), bean2.getName());

		// exact size
		if (CHECK_EXACT_OUTPUT_SIZE) {
			Assert.assertEquals(bean1.getContentLength(), bean2.getContentLength());
		}

		// zero size not allowed if source non-zero
		if (bean1.getContentLength() > 0) {
			Assert.assertTrue(bean2.getContentLength() > 0, "zero size dataset");
		}
	}

	
	public static void main(String[] args) throws Exception {

		boolean testOK = false;
		try {

			String configURL = null;
			String username = null;
			String password = null;
			String sessionsDir = DEFAULT_SESSIONS_DIR;
			switch(args.length) {
			case 3:
				configURL = args[0];
				username = args[1];
				password = args[2];
				DirectoryLayout.initialiseClientLayout(configURL);
				break;
			case 4:
				configURL = args[0];
				username = args[1];
				password = args[2];
				DirectoryLayout.initialiseClientLayout(configURL);
				sessionsDir = args[3];
				break;

			default:
				System.out.println("Usage: " + SessionReplayTest.class.getSimpleName() + " <config-url username password> <sessions dir>\n" +
				"If there are no arguments, config file is used.");
				updateFlagFileAndExit(false);
			}


			// initialize
			SessionReplayTest test = null;
			try {
				test = new SessionReplayTest(username, password, configURL);
			} catch (Exception e) {
				e.printStackTrace();
				System.out.println("TOOL TESTS INIT ERROR");
				updateFlagFileAndExit(false);
			}

			// run
			try {
				test.setUp();
				testOK = test.testSessions(sessionsDir);
			} catch (Exception e) {
				e.printStackTrace();
				System.out.println("TOOL TEST ERROR");
				test.tearDown();
				updateFlagFileAndExit(false);
			} finally {
				test.tearDown();
			}

			if (testOK) {
				System.out.println("TOOL TESTS OK");
				updateFlagFileAndExit(true);
			} else {
				System.out.println("TOOL TESTS FAILED");
				updateFlagFileAndExit(false);
			}
		} 
		
		// Should only get here if something weird happens
		finally {
			updateFlagFileAndExit(false);
		}
	}
	
	private void createReports() throws IOException {
		
		
		// count unique tools and sessions
		HashMap<String, Integer> uniqueTools = new HashMap<String,Integer>(); 
		Set<File> uniqueSessions = new HashSet<File>();
		HashMap<String, List<File>> toolToSessionsMap = new HashMap<String, List<File>>();
		for (ToolTestResult toolTestResult : toolTestResults) {
			if (!uniqueTools.containsKey(toolTestResult.getOperation().getID())) {
				uniqueTools.put(toolTestResult.getOperation().getID(), 1);
			} else {
				uniqueTools.put(toolTestResult.getOperation().getID(), uniqueTools.get(toolTestResult.getOperation().getID()) + 1);
			}

			if (!toolToSessionsMap.containsKey(toolTestResult.getOperation().getID())) {
				List<File> sessionsList = new LinkedList<File>();
				sessionsList.add(toolTestResult.getSession());
				toolToSessionsMap.put(toolTestResult.getOperation().getID(), sessionsList);
			} else {
				List<File> sessionsList = toolToSessionsMap.get(toolTestResult.getOperation().getID());
				if (!sessionsList.contains(toolTestResult.getSession())) {
					sessionsList.add(toolTestResult.getSession());
				}
			}
			
			if (!uniqueSessions.contains(toolTestResult.getSession())) {
				uniqueSessions.add(toolTestResult.getSession());
			}
		}

		// sort tools by test count
		TreeMap<Integer,String> testsPerTool = new TreeMap<Integer,String>();
		for (Entry<String, Integer> entry : uniqueTools.entrySet()) {
			System.out.println("putting " + entry.getValue() + " " + entry.getKey());
			testsPerTool.put(new Integer(entry.getValue()), entry.getKey());
		}
		
		// failed and successful tasks
		List<ToolTestResult> failedTasks = new LinkedList<ToolTestResult>();
		List<ToolTestResult> successTasks = new LinkedList<ToolTestResult>();
		for (ToolTestResult toolTestResult : toolTestResults) {
			if (TestResult.FAIL.equals(toolTestResult.getTestResult())) {
				failedTasks.add(toolTestResult);
			} else {
				successTasks.add(toolTestResult);
			}
		}
		
		//

		File htmlFile = new File("index.html");
		FileWriter writer = new FileWriter(htmlFile);
		writer.write("<html><body>");
		writer.write("<h2>Tool tests</h2>");

		writer.write("<h3>Summary</h3>");
		writer.write("<p>Tasks: " + failedTasks.size() + " failed " + successTasks.size() + " successful "+ toolTestResults.size() + " total</p>");

		writer.write("<p>Unique tools tested: " + uniqueTools.size() + "/" + getTotalNumberOfTools() + "</p>");
		writer.write("<p>Number of sessions: " + uniqueSessions.size() + "</p>");

		
		// Failed tests
		writer.write("<h3>Failed tools</h3>");
		writer.write("<table><tr>" + 
				"<th>Tool</th>" + 
				"<th>Session</th>" +
				"<th>Test error message</th>" + 
				"<th>Task end state</th>" + 
				"<th>Task error message</th>" +
				"<th>Task screen output</th>" + 
				"</tr>");
		for (ToolTestResult toolTestResult : failedTasks) {
			writer.write("<tr>" +
					"<td>" + toolTestResult.getOperation().getDefinition().getFullName() + "</td>" +
					"<td>" + toolTestResult.getSession().getName() + "</td>" + // TODO add link
					"<td>" + toolTestResult.getTestErrorMessage() + "</td>" +
					"<td>" + toolTestResult.getTask().getState() + "</td>" +
					"<td>" + toolTestResult.getTask().getErrorMessage() + "</td>" +
//					"<td>" + toolTestResult.getTask().getScreenOutput() + "</td>" +
					"<td>" + "link to output" + "</td>" +
			    	"</tr>");
		}
		writer.write("</table>");

		// Successful tools
		writer.write("<h3>Successful tools</h3>");
		writer.write("<table><tr>" + 
				"<th>Tool</th>" + 
				"<th>Session</th>" +
				"<th>Test error message</th>" + 
				"<th>Task end state</th>" + 
				"<th>Task error message</th>" +
				"<th>Task screen output</th>" + 
				"</tr>");
		for (ToolTestResult toolTestResult : successTasks) {
			writer.write("<tr>" +
					"<td>" + toolTestResult.getOperation().getDefinition().getFullName() + "</td>" +
					"<td>" + toolTestResult.getSession().getName() + "</td>" + // TODO add link
					"<td>" + toolTestResult.getTestErrorMessage() + "</td>" +
					"<td>" + toolTestResult.getTask().getState() + "</td>" +
					"<td>" + toolTestResult.getTask().getErrorMessage() + "</td>" +
					"<td>" + "link to output" + "</td>" +
//					"<td>" + toolTestResult.getTask().getScreenOutput() + "</td>" +
			    	"</tr>");
		}
		writer.write("</table>");
		
		
		writer.write("<h3>Coverage</h3>");
		writer.write("<table><tr>" + 
				"<th>Tool</th>" + 
				"<th>Test count</th>" +
				"<th>Sessions</th>" + 
				"</tr>");
		for (Entry<Integer, String> entry : testsPerTool.entrySet()) {
			String sessionsString = "";
			for (File session : toolToSessionsMap.get(entry.getValue())) {
				sessionsString += session.getName() + " ";
			}
			sessionsString.trim();
			
			writer.write("<tr>" +
					"<td>" + this.getOperationDefinition(entry.getValue(), toolModules).getFullName()  + "</td>" +
					"<td>" + entry.getKey() + "</td>" +
					"<td>" + sessionsString + "</td>" +
			    	"</tr>");
		}
		writer.write("</table>");

		
		writer.write("</body></html>");
		writer.flush();
		writer.close();
	}
	
	private static void updateFlagFileAndExit(boolean testOK) throws IOException {
		File flagFile = new File(FLAG_FILE);
		if (testOK) {
			if (!flagFile.exists()) {
				flagFile.createNewFile();
			} else {
				flagFile.setLastModified(System.currentTimeMillis());
			}
		} else {
			flagFile.delete();
		}
		
		if (testOK) {
			System.exit(0);
		} else {
			System.exit(1);
		}
	}
	
	
	private enum TestResult {
		OK,FAIL;
	}
	
	private class ToolTestResult {
		private TestResult testResult;
		private Operation operation;
		private File sessionFile;
		private Task task;
		private String testErrorMessage;
		
		public ToolTestResult(TestResult testResult, Operation operation, File sessionFile, Task task) {
			this(testResult, operation, sessionFile, task, null);
		}
		
		public ToolTestResult(TestResult testResult, Operation operation, File sessionFile, Task task, String testErrorMessage) {
			this.testResult = testResult;
			this.operation = operation;
			this.sessionFile = sessionFile;
			this.task = task;
			this.testErrorMessage = testErrorMessage;
		}

		public TestResult getTestResult() {
			return testResult;
		}
		public Operation getOperation() {
			return operation;
		}
		public File getSession() {
			return sessionFile;
		}
		public Task getTask() {
			return task;
		}
		
		public String getTestErrorMessage() {
			return testErrorMessage;
		}
		
	}
	
	private int getTotalNumberOfTools() {
		Set<String> uniqueIDs = new HashSet<String>();
		for (ToolModule toolModule : toolModules) {
			for (ToolCategory toolCategory : toolModule.getVisibleCategories()) {
				for (OperationDefinition od : toolCategory.getToolList()) {
					if (!uniqueIDs.contains(od.getID())) {
						uniqueIDs.add(od.getID());
					}
				}
			}
		}
		return uniqueIDs.size();
	}
	
	
	public static class SessionLoadingSkeletonApplication extends ClientApplication {

		private SessionReplayTest parent;

		public SessionLoadingSkeletonApplication(SessionReplayTest parent, LinkedList<ToolModule> toolModules) {
			this.parent = parent;
			this.toolModules = toolModules;
			logger = org.apache.log4j.Logger.getLogger(SessionLoadingSkeletonApplication.class);
		}
		
		@Override
		public OperationDefinition getOperationDefinition(String toolId) {
			return parent.getOperationDefinition(toolId, toolModules);
		};
		
		@Override
		public void checkFreeMemory() {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void createLink(DataBean source, DataBean target, Link type) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void deleteDatas(DataItem... datas) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void flipTaskListVisibility(boolean closeIfVisible) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public DataManager getDataManager() {
			return manager;
//			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public Icon getIconFor(DataItem data) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public VisualisationFrameManager getVisualisationFrameManager() {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void importGroup(Collection<ImportItem> datas, String folderName) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		protected void initialiseGUI() throws MicroarrayException, IOException {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public DataFolder initializeFolderForImport(String folderName) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void loadSession() {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void loadSessionFrom(URL url) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void loadSessionFrom(File file) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public File openWorkflow() {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void removeLink(DataBean source, DataBean target, Link type) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void reportException(Exception e) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void reportInitialisation(String report, boolean newline) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void reportTaskError(Task job) throws MicroarrayException {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void restoreSessionFrom(File file) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void runBlockingTask(String taskName, Runnable runnable) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void runWorkflow(URL workflowScript) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void runWorkflow(URL workflowScript, AtEndListener atEndListener) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void saveSession() {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public File saveWorkflow() {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void setMaximisedVisualisationMode(boolean maximisedVisualisationMode) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void showDetailsFor(DataBean data) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void showDialog(String title, String message, String details, Severity severity, boolean modal) {
			showDialog(title, message, details, severity, modal, DetailsVisibility.DETAILS_VISIBLE, null);
		}

		@Override
		public void showDialog(String title, String message, String details, Severity severity, boolean modal, DetailsVisibility detailsVisibility, PluginButton button) {
			showDialog(title, message, details, severity, modal, detailsVisibility, button, false);
		}

		@Override
		public void showDialog(String title, String message, String details, Severity severity, boolean modal, DetailsVisibility detailsVisibility, PluginButton button, boolean feedBackEnabled) {
			System.out.println(title + "\n" + message + "\n" + details);
		}

		@Override
		public void showHistoryScreenFor(DataBean data) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void showImportToolFor(File file, String destinationFolder, boolean skipActionChooser) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void showPopupMenuFor(MouseEvent e, DataItem data) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void showPopupMenuFor(MouseEvent e, List<DataItem> datas) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void showSourceFor(String operationName) throws TaskException {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		protected void taskCountChanged(int newTaskCount, boolean attractAttention) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void viewHelp(String id) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void viewHelpFor(OperationDefinition operationDefinition) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void visualiseWithBestMethod(FrameType target) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}
		
	}

}
