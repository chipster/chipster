package fi.csc.microarray.analyser;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URL;
import java.text.ParseException;
import java.util.Date;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.UUID;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import org.junit.Assert;
import org.springframework.validation.Errors;

import com.beust.jcommander.JCommander;

import fi.csc.microarray.analyser.AnalysisTestBase.JobResultListener;
import fi.csc.microarray.analyser.SessionReplayTestParameters.CommandGenerate;
import fi.csc.microarray.analyser.SessionReplayTestParameters.CommandRun;
import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.RemoteServiceAccessor;
import fi.csc.microarray.client.ServiceAccessor;
import fi.csc.microarray.client.Session;
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
import fi.csc.microarray.client.session.SessionManager;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.tasks.Task.State;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.DataNotAvailableHandling;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.DataManager.ContentLocation;
import fi.csc.microarray.databeans.DataManager.StorageMethod;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.MessagingTestBase;
import fi.csc.microarray.module.ModuleManager;
import fi.csc.microarray.module.chipster.MicroarrayModule;
import fi.csc.microarray.util.Exceptions;
import fi.csc.microarray.util.IOUtils;


public class SessionReplayTest extends MessagingTestBase {

	public static final String SESSIONS_DIR_NAME = "sessions";
	public static final String SCREEN_OUTPUTS_DIR_NAME = "screen-outputs";
	public static final String STACKTRACES_DIR_NAME = "stacktraces";
	
	private static final String FLAG_FILE = "test-ok";
	private static final String DEFAULT_WEB_DIR = "web";
	public static final String SCREEN_OUTPUT_POSTFIX = "-output.txt";
	public static final String STACKTRACE_POSTFIX = "-stacktrace.txt";
	
	private static final long TOOL_TEST_TIMEOUT = 3;
	private static final TimeUnit TOOL_TEST_TIMEOUT_UNIT = TimeUnit.HOURS;
	
	private static final boolean FAIL_ON_OUTPUT_SIZE_MISMATCH = false;
	private static final boolean CHECK_CONTENTS = true;
	
	private File sessionsDir;
	static private File webDir; // needed in early fail
	private File screenOutputsDir;
	private File stacktracesDir;
	
	
	private TaskExecutor executor;
	private DataManager manager, sourceManager;
	private ServiceAccessor serviceAccessor;
	LinkedList<ToolModule> toolModules;
	
	public SessionReplayTest(String username, String password, String configURL) {
		super(username, password, configURL);

		screenOutputsDir = new File(webDir, SCREEN_OUTPUTS_DIR_NAME);
		screenOutputsDir.mkdirs();
		
		// delete existing screen outputs
		if (screenOutputsDir.exists()) {
			for (File f : screenOutputsDir.listFiles()) {
				f.delete();
			}
		}

		stacktracesDir = new File(webDir, STACKTRACES_DIR_NAME);
		stacktracesDir.mkdirs();
		
		// delete existing
		if (stacktracesDir.exists()) {
			for (File f : stacktracesDir.listFiles()) {
				f.delete();
			}
		}
	}

	public boolean testSessions(String sessionsDirName) throws Exception {

		this.sessionsDir = new File(sessionsDirName);
		ToolTestSummary summary = new ToolTestSummary();
		summary.setStartTime(new Date());
		
		// Set up modules
		ModuleManager moduleManager = new ModuleManager("fi.csc.microarray.module.chipster.MicroarrayModule");
		Session.getSession().setModuleManager(moduleManager);
		try {
			// Set up main (target) system
			manager = new DataManager();
			moduleManager.plugAll(manager, null);
			executor = new TaskExecutor(super.endpoint, manager);
			toolModules = new LinkedList<ToolModule>();
			serviceAccessor = new RemoteServiceAccessor();
			serviceAccessor.initialise(manager, this.authenticationListener);
			serviceAccessor.fetchDescriptions(new MicroarrayModule());
			Session.getSession().setServiceAccessor(serviceAccessor);
			toolModules.addAll(serviceAccessor.getModules());
			Session.getSession().setClientApplication(new SessionLoadingSkeletonApplication(this, toolModules, manager));

			// Set up source system
			sourceManager = new DataManager();
			moduleManager.plugAll(sourceManager, null);

			// Run all sessions
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
					testSession(testSession, summary);
				} catch (Throwable e) {
					e.printStackTrace();
					File f = writeThrowableToDisk(e);
					summary.getSessionsWithErrors().put(testSession, f.getName());
				}
			}
		} finally {
			if (serviceAccessor != null) {
				serviceAccessor.close();
			}
		}

		// running tests done

		
		// create result files
		addAllToolsToSummary(summary);
		summary.calculateStats();
		summary.setEndTime(new Date());
		summary.writeToFiles(webDir);
		
		return getOverallResult(summary);
		
	}

	private static boolean getOverallResult(ToolTestSummary summary) {
		
		// Get overall result
		boolean combinedToolsResult = true;
		for (ToolTestResult toolTestResult : summary.getToolTestResults()) {
			if (toolTestResult.getTestResult().equals(TestResult.FAIL)) {
				combinedToolsResult = false;
			}
		}
		if (!combinedToolsResult) {
			System.out.println("failed tools, failing");
			return false;
		}
		
		// Fail if no results at all
		if (summary.getToolTestResults().size() == 0) {
			System.out.println("zero results, failing");
			return false;
		}
		
		// Fail if missing tools
		if (summary.getSessionsWithMissingTools().size() > 0) {
			System.out.println("missing tools, failing");
			return false;
		}

		// Fail if sessions with errors
		if (summary.getSessionsWithErrors().size() > 0) {
			System.out.println("sessions with errors, failing");
			return false;
		}
		
		return true;

	}
	
	private void testSession(File session, ToolTestSummary summary) throws Exception {

		Map<DataBean, DataBean> sourceDataBeanToTargetDataBean = new HashMap<DataBean, DataBean>();

		// Load session
		// some loading code uses Session.getSession().getDataManagers
		Session.getSession().setDataManager(sourceManager);
		SessionManager sourceSessionManager = new SessionManager(sourceManager, serviceAccessor.getTaskExecutor(), serviceAccessor.getFileBrokerClient(), null);
		sourceSessionManager.loadLocalSession(session, false);
		Session.getSession().setDataManager(manager);
		
		// Pick import operations and copy imported data beans to target manager 
		// Also map OperationRecords to outputs TODO check that order is right, might need to traverse links
		LinkedList<OperationRecord> rootLevelOperationRecords = new LinkedList<OperationRecord>();
		Map<OperationRecord, List<DataBean>> outputMap = new HashMap<OperationRecord, List<DataBean>>();
		LinkedList<DataBean> dataBeansInSourceManagerWhichWereCopied = new LinkedList<DataBean>();
		for (DataBean dataBean : sourceManager.databeans()) {
			OperationRecord operationRecord = dataBean.getOperationRecord();

			// pick import operations, local operations 
			// and those which did have inputs when run, but don't have parents now (i.e. inputs have
			// been deleted for example to save space)
			if (OperationDefinition.IMPORT_DEFINITION_ID.equals(operationRecord.getNameID().getID()) ||
					"LocalNGSPreprocess.java".equals(operationRecord.getNameID().getID()) ||
					(dataBean.getLinkTargets(Link.derivationalTypes()).size() == 0 &&
					operationRecord.getInputRecords().size() > 0)) {
				
				// load imported databean, add mapping
				DataBean dataBeanCopy = manager.createDataBean(dataBean.getName());
				URL urlInSessionZip = null;
				for (ContentLocation contentLocation : sourceManager.getContentLocationsForDataBeanSaving(dataBean)) {
					if (contentLocation.getMethod() == StorageMethod.LOCAL_SESSION_ZIP) {
						urlInSessionZip = contentLocation.getUrl(); 
					}
				}
				
				if (urlInSessionZip == null) {
					throw new IllegalArgumentException("session file " + session.getName() + " must contain all data files (missing " + dataBean.getName() + ")");
				}
				URL url = new URL(session.toURI().toURL(), "#" + urlInSessionZip.getRef());
				manager.addContentLocationForDataBean(dataBeanCopy, StorageMethod.LOCAL_SESSION_ZIP, url);

				sourceDataBeanToTargetDataBean.put(dataBean, dataBeanCopy);
				dataBeansInSourceManagerWhichWereCopied.add(dataBean);
				
				// avoid NPE 
				dataBeanCopy.setOperationRecord(operationRecord);
				
				// TODO what if not in the root folder in the source manager
				manager.connectChild(dataBeanCopy, manager.getRootFolder());
				rootLevelOperationRecords.add(operationRecord);
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

		
		// Add links between copied (imported, without parents) data beans
		for (DataBean originalBean : dataBeansInSourceManagerWhichWereCopied) {
			for (Link linkType : Link.values()) {
				for (DataBean originalLinkSourceBean : originalBean.getLinkSources(linkType)) {
					if (dataBeansInSourceManagerWhichWereCopied.contains(originalLinkSourceBean)) {
						sourceDataBeanToTargetDataBean.get(originalLinkSourceBean).addLink(linkType, sourceDataBeanToTargetDataBean.get(originalBean));
					}
				}
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
			if (rootLevelOperationRecords.contains(operationRecord)) {
				System.out.println("skipping root level operation " + operationRecord.getFullName());
				continue;
			}

			System.out.println("setting up " + operationRecord.getFullName());

			// Get inputs
			LinkedList <DataBean> inputBeans = new LinkedList<DataBean>();
			for (InputRecord inputRecord : operationRecord.getInputRecords()) {
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
			OperationDefinition operationDefinition = getOperationDefinition(operationRecord.getNameID().getID(), toolModules);
			if (operationDefinition == null) {
				System.out.println("Missing tool: " + operationRecord.getNameID().getID() + "  in session: " + session.getName() + " skipping rest of the session");
				summary.getSessionsWithMissingTools().put(session, operationRecord.getNameID().getID());
				return;
			}
			Operation operation = new Operation(operationDefinition, inputBeans.toArray(new DataBean[] {}));

			// Parameters, copy paste from workflows
			for (ParameterRecord parameterRecord : operationRecord.getParameters()) {
				if (parameterRecord.getValue() != null && !parameterRecord.getValue().equals("")) {	
					Parameter definitionParameter = operation.getDefinition().getParameter(parameterRecord.getNameID().getID()); 
					if (definitionParameter != null) {
						Parameter parameter = (Parameter)definitionParameter.clone();
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

			Task task = executor.createNewTask(new OperationRecord(operation), operation.getDefinition().isLocal());
			
			// Execute the task
			System.out.println(new Date() + " running " + operation.getDefinition().getFullName());
			CountDownLatch latch = new CountDownLatch(1);
			task.addTaskEventListener(new JobResultListener(latch));
			executor.startExecuting(task);
			latch.await(TOOL_TEST_TIMEOUT, TOOL_TEST_TIMEOUT_UNIT);
			
			// write screen output to disk
			writeScreenOutputToDisk(task);
			
			// Task failed
			if (!task.getState().equals(State.COMPLETED)) {
				
				// try to cancel if still running
				if (task.getState().equals(State.RUNNING)) {
					executor.kill(task);
					summary.getToolTestResults().add(new ToolTestResult(TestResult.FAIL, session, task, "task did not finish before test timeout " +TOOL_TEST_TIMEOUT + " " + TOOL_TEST_TIMEOUT_UNIT.toString()));
					return;
				}
				summary.getToolTestResults().add(new ToolTestResult(TestResult.FAIL, session, task, "task was not completed"));
				return;
			}
			
			// Check results
			// Target manager needs to be available through session for some of these to work
			Session.getSession().setDataManager(manager);
			try {
				
				// Link result beans, add to folders etc			
				Session.getSession().getApplication().onFinishedTask(task, operation.getResultListener(), task.getState(), true);

				// Check that number of results and result names match
				Iterator<DataBean> targetIterator = task.getOutputs().iterator();
				for (DataBean sourceBean : outputMap.get(operationRecord)) {
					if (targetIterator.hasNext()) {
						DataBean targetBean = targetIterator.next();

						// check content types
						if (!sourceBean.getContentType().getType().equals(targetBean.getContentType().getType())) {
							summary.getToolTestResults().add(new ToolTestResult(TestResult.FAIL, session, task, "Mismatch in result content types, "
									+ sourceBean.getName() + ": " + sourceBean.getContentType().getType() + ", " + targetBean.getName() + ": " + targetBean.getContentType().getType()));
							return;
						}
						
//						// check names
//						if (!sourceBean.getName().equals(targetBean.getName())) {
//							toolTestResults.add(new ToolTestResult(TestResult.FAIL, operation, session, task, "mismatch in result dataset names, "
//									+ "expecting: " + sourceBean.getName() + " got: " + targetBean.getName()));
//							return;
//						}
					} 
					// Not enough results
					else {
						summary.getToolTestResults().add(new ToolTestResult(TestResult.FAIL, session, task, "not enough result datasets"));
						return;
					}
				}
				if (targetIterator.hasNext()) {
					summary.getToolTestResults().add(new ToolTestResult(TestResult.FAIL, session, task, "too many result datasets"));
					return;
				}

				// Find and replace metadata 
				targetIterator = task.getOutputs().iterator();
				for (DataBean sourceBean : outputMap.get(operationRecord)) {
					DataBean targetBean = targetIterator.next();

					// replace metadata contents from the source session
					if (Session.getSession().getPrimaryModule().isMetadata(targetBean)) {
						System.out.println("copying metadata for: " + targetBean.getName());
						OutputStream metadataOut = manager.getContentOutputStreamAndLockDataBean(targetBean);
						InputStream sourceIn = null;
						try {
							sourceIn = sourceManager.getContentStream(sourceBean, DataNotAvailableHandling.EXCEPTION_ON_NA);
							IOUtils.copy(sourceIn, metadataOut);
						} finally {
							IOUtils.closeIfPossible(sourceIn);
							manager.closeContentOutputStreamAndUnlockDataBean(targetBean, metadataOut);
						}
					}
				}

				// Compare data beans, add source bean -> target bean mapping
				List<String> outputsWithMisMatchingSizes = new LinkedList<String>();
				List<String> outputsWithMisMatchingContents = new LinkedList<String>();

				targetIterator = task.getOutputs().iterator();
				for (DataBean sourceBean : outputMap.get(operationRecord)) {
					DataBean targetBean = targetIterator.next();
					try {
						compareDataBeans(sourceBean, targetBean);
					} catch (Throwable t) {
						summary.getToolTestResults().add(new ToolTestResult(TestResult.FAIL, session, task, t.getMessage()));
						return;
					}
					
					// Add source bean -> target bean mapping, needed for further operations
					sourceDataBeanToTargetDataBean.put(sourceBean, targetBean);
					
					// Collect size matches
					if (Session.getSession().getDataManager().getContentLength(sourceBean) != Session.getSession().getDataManager().getContentLength(targetBean)) {
						if (sourceBean.getName().equals(targetBean.getName())) {
							outputsWithMisMatchingSizes.add(sourceBean.getName());
						} else {
							outputsWithMisMatchingSizes.add(sourceBean.getName() + " | " + targetBean.getName());
						}
					}

					// Collect content matches
					if (CHECK_CONTENTS) {
						InputStream sourceIn = null, targetIn = null;
						try {
							sourceIn = sourceManager.getContentStream(sourceBean, DataNotAvailableHandling.EXCEPTION_ON_NA);
							targetIn = manager.getContentStream(targetBean, DataNotAvailableHandling.EXCEPTION_ON_NA);
							if (!IOUtils.contentEquals(sourceIn, targetIn)) {

								if (sourceBean.getName().equals(targetBean.getName())) {
									outputsWithMisMatchingContents.add(sourceBean.getName());
								} else {
									outputsWithMisMatchingContents.add(sourceBean.getName() + " | " + targetBean.getName());
								}
							}
						} finally {
							IOUtils.closeIfPossible(sourceIn);
							IOUtils.closeIfPossible(targetIn);
						}
					}
				}

				// Add result
				ToolTestResult toolTestResult = new ToolTestResult(TestResult.OK, session, task, null);
				toolTestResult.setOutputsWithMisMatchingSizes(outputsWithMisMatchingSizes);
				toolTestResult.setOutputsWithMisMatchingContents(outputsWithMisMatchingContents);
				summary.getToolTestResults().add(toolTestResult);

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

//		// Check name, may be already checked though
//		Assert.assertEquals(bean1.getName(), bean2.getName());

		// check content types
		Assert.assertEquals(bean1.getContentType().getType(), bean2.getContentType().getType());
		
		// exact size
		if (FAIL_ON_OUTPUT_SIZE_MISMATCH) {
			Assert.assertEquals("comparing output size", Session.getSession().getDataManager().getContentLength(bean1), Session.getSession().getDataManager().getContentLength(bean2));
		}

		// zero size not allowed if source non-zero
		if (Session.getSession().getDataManager().getContentLength(bean1) > 0) {
			Assert.assertTrue("zero size dataset", Session.getSession().getDataManager().getContentLength(bean2) > 0);
		}
	}

	
	public static void main(String[] args) throws Exception {

		// parse command line parameters
		SessionReplayTestParameters parameters = new SessionReplayTestParameters();
		JCommander jc = new JCommander(parameters);
		CommandRun run = new CommandRun();
		jc.addCommand("run", run);
		CommandGenerate generate = new CommandGenerate();
		jc.addCommand("generate", generate);
		try {
			jc.parse(args);
		} catch (Exception e) {
			System.out.println(e.getMessage());
			jc.usage();
			initWebDir(null);
			System.out.println("TOOL TEST EARLY ERROR");
			updateFlagFileAndExit(false);
		}
		
		try { 
			if ("run".equals(jc.getParsedCommand())) {
				System.out.println(new Date() + " running session replay test");
				initWebDir(run.output);
				runTestSessions(run.config, run.username, run.password, run.sessions);
				
				
			} else if ("generate".equals(jc.getParsedCommand())) {
				initWebDir(generate.output);
				File resultsDir = new File(generate.generateFrom);
				readResultsFromFiles(resultsDir);
			
			} else {
				jc.usage();
				initWebDir(null);
				updateFlagFileAndExit(false);
			}
		} catch (Exception e) {
			System.out.println(Exceptions.getStackTrace(e));
			System.out.println("TOOL TEST EXCEPTION ERROR");
			updateFlagFileAndExit(false);
		}

		System.out.println("TOOL TEST UNEXCPECTED ERROR");
		updateFlagFileAndExit(false);
			
	}
	
	private static void initWebDir(String dirName) {
		if (dirName == null || dirName.isEmpty()) {
			webDir = new File(DEFAULT_WEB_DIR); 
		} else {
			webDir = new File(dirName);
		}
		webDir.mkdirs();
	}
	

	private static void runTestSessions(String config, String username, String password, String sessionsDirName) throws Exception {

		SessionReplayTest test = null;
		try {
			DirectoryLayout.initialiseClientLayout(config);
			test = new SessionReplayTest(username, password, config);
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("TOOL TESTS INIT ERROR");
			updateFlagFileAndExit(false);
		}

		try {
			test.setUp();
			// run tests
			boolean testOK = test.testSessions(sessionsDirName);
			test.tearDown();
			if (testOK) {
				System.out.println("TOOL TESTS OK");
				updateFlagFileAndExit(true);
			} else {
				System.out.println("TOOL TESTS FAILED");
				updateFlagFileAndExit(false);
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("TOOL TEST ERROR");
			test.tearDown();
			updateFlagFileAndExit(false);
		} finally {
			// never get here?
			test.tearDown();
		}
	}
	
	
	private static void readResultsFromFiles(File inputDir) throws IOException, ParseException {
		ToolTestSummary summary = new ToolTestSummary();
		summary.readFromFiles(inputDir);
		summary.calculateStats();
		summary.writeToFiles(webDir);
		boolean summaryResult = getOverallResult(summary);
		System.out.println("summary result is : " + summaryResult);
		updateFlagFileAndExit(summaryResult);
	}
	
	private void writeScreenOutputToDisk(Task task) {
		if (task.getScreenOutput() != null) {
			File outputFile = new File(screenOutputsDir, task.getId() + SCREEN_OUTPUT_POSTFIX);
			try {
				IOUtils.copy(new ByteArrayInputStream(task.getScreenOutput().getBytes()), outputFile);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}	
	
	private File writeThrowableToDisk(Throwable t) {
		File outputFile = new File(stacktracesDir, UUID.randomUUID().toString() + STACKTRACE_POSTFIX );
		try {
			IOUtils.copy(new ByteArrayInputStream(Exceptions.getStackTrace(t).getBytes()), outputFile);
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
		return outputFile;
	}

	
	private static void updateFlagFileAndExit(boolean testOK) throws IOException {
		File flagFile = new File(webDir, FLAG_FILE);
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
	
	
	enum TestResult {
		OK,FAIL;
	}
	
	private void addAllToolsToSummary(ToolTestSummary summary) {
		for (ToolModule toolModule : toolModules) {
			for (ToolCategory toolCategory : toolModule.getVisibleCategories()) {
				for (OperationDefinition od : toolCategory.getToolList()) {
					summary.getAllTools().put(od.getID(), od.getFullName());
				}
			}
		}
	}

	
	
	public static class SessionLoadingSkeletonApplication extends ClientApplication {

		private SessionReplayTest parent;

		public SessionLoadingSkeletonApplication(SessionReplayTest parent, LinkedList<ToolModule> toolModules, DataManager manager) {
			this.parent = parent;
			this.toolModules = toolModules;
			this.manager = manager;
			logger = org.apache.log4j.Logger.getLogger(SessionLoadingSkeletonApplication.class);
		}
		
		@Override
		public OperationDefinition getOperationDefinition(String toolId) {
			return parent.getOperationDefinition(toolId, toolModules);
		};

		@Override
		public void initialiseGUIThreadSafely(File session) throws MicroarrayException, IOException {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void reportException(Exception e) {
			System.err.println(this.getClass().getSimpleName() + " reporting an exception:");
			e.printStackTrace();
			logger.error(e);
		}

		@Override
		public void reportInitialisationThreadSafely(String report, boolean newline) {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void reportTaskError(Task job) throws MicroarrayException {
			throw new UnsupportedOperationException("not supported by skeleton app");
		}

		@Override
		public void runBlockingTask(String taskName, Runnable runnable) {
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
		public void reportExceptionThreadSafely(Exception e) {
			System.err.println(this.getClass().getSimpleName() + " reporting an exception:");
			e.printStackTrace();
			logger.error(e);
		}
	}


	
}
