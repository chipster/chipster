package fi.csc.microarray.analyser;

import java.awt.event.MouseEvent;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URL;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import javax.jms.JMSException;
import javax.swing.Icon;

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
import fi.csc.microarray.client.operation.ToolModule;
import fi.csc.microarray.client.operation.OperationRecord.InputRecord;
import fi.csc.microarray.client.operation.OperationRecord.ParameterRecord;
import fi.csc.microarray.client.operation.parameter.DataSelectionParameter;
import fi.csc.microarray.client.operation.parameter.Parameter;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.tasks.TaskException;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.client.tasks.Task.State;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.MessagingTestBase;
import fi.csc.microarray.module.ModuleManager;
import fi.csc.microarray.module.chipster.MicroarrayModule;
import fi.csc.microarray.util.IOUtils;

public class SessionReplayTest extends MessagingTestBase {

	private static final boolean CHECK_EXACT_OUTPUT_SIZE = true;
	
	private File[] testSessions = new File[] {
			new File("/home/hupponen/ngs-test.zip"),
			new File("/home/hupponen/ngs-test-many.zip"),
			new File("/home/hupponen/ngs-test-with-changed-parameter.zip"),
			new File("/home/hupponen/test-example.zip"),

	};
	
	public SessionReplayTest(String username, String password, String configURL) {
		super(username, password, configURL);
	}

	public boolean test() throws Exception {
		LinkedHashMap<File, Boolean> sessionTestResults = new LinkedHashMap<File, Boolean>();
		boolean combinedResult = true;
		for (File testSession : testSessions) {

			boolean result = false;
			try {
				result = testSession(testSession);
			} catch (Throwable e) {
				result = false;
			}
			sessionTestResults.put(testSession, result);
			if (!result) {
				combinedResult = false;
			}
		
		}

		System.out.println("-----------------------------");
		for (File testSession : sessionTestResults.keySet()) {
			String result = sessionTestResults.get(testSession) ? "OK" : "FAILED";
			System.out.println(testSession.getName() + " " +  result);
		}
		System.out.println("-----------------------------");
		
		return combinedResult;
	}

	private boolean testSession(File session) throws InstantiationException,
			IllegalAccessException, ClassNotFoundException, IOException,
			JMSException, Exception, MicroarrayException, TaskException,
			InterruptedException {
		// Set up modules
		ModuleManager moduleManager = new ModuleManager("fi.csc.microarray.module.chipster.MicroarrayModule");
		Session.getSession().setModuleManager(moduleManager);
		
		// Set up main (target) system
		DataManager manager = new DataManager();
		moduleManager.plugAll(manager, null);
		TaskExecutor executor = new TaskExecutor(super.endpoint, manager);
		LinkedList<ToolModule> toolModules = new LinkedList<ToolModule>();
		ServiceAccessor serviceAccessor = new RemoteServiceAccessor();
		serviceAccessor.initialise(manager, this.authenticationListener);
		serviceAccessor.fetchDescriptions(new MicroarrayModule());
		toolModules.addAll(serviceAccessor.getModules());
		Session.getSession().setClientApplication(new SessionLoadingSkeletonApplication(this, toolModules));
		
		// Set up source system
		DataManager sourceManager = new DataManager();
		moduleManager.plugAll(sourceManager, null);
		sourceManager.loadSession(session, false);
		
		Map<DataBean, DataBean> sourceDataBeanToTargetDataBean = new HashMap<DataBean, DataBean>();
		
		
		// Pick import operations and copy imported data beans to target manager 
		// Also map OperationRecords to outputs TODO check that order is right, might need to traverse links
		LinkedList<OperationRecord> importOperationRecords = new LinkedList<OperationRecord>();
		Map<OperationRecord, List<DataBean>> outputMap = new HashMap<OperationRecord, List<DataBean>>();
		for (DataBean dataBean : sourceManager.databeans()) {
			OperationRecord operationRecord = dataBean.getOperationRecord();

			// pick import operations
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
		
		System.out.println("source session databeans: " + sourceManager.databeans().size());
		System.out.println("source session operation records: " + operationRecords.size());
		
		

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
					System.out.println("not enough inputs for " + operationRecord.getFullName() + ", skipping");
					// TODO fail the test
					continue;
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

						// set it
						//System.out.println("setting parameter " + parameter.getID() + ": " + parameter.getValueAsString());
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
			latch.await(1000*60*10, TimeUnit.MILLISECONDS); // FIXME remove hard coding
			if (!task.getState().equals(State.COMPLETED)) {
				return false;
			}
			
			// Link data beans, link metadata etc
			// Target manager needs to be available throught session for some of these to work
			Session.getSession().setDataManager(manager);
			try {
				Session.getSession().getApplication().onFinishedTask(task, operation);

				// Add source bean -> target bean mapping
				// There might be more target outputs than source outputs, it's fine
				Iterator<DataBean> targetIterator = task.outputs().iterator();
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

					Assert.assertEquals(sourceBean.getName(), targetBean.getName());
					sourceDataBeanToTargetDataBean.put(sourceBean, targetBean);
				}
			} finally {
				// Set session data manager back to null to avoid problems
				Session.getSession().setDataManager(null);
			}

			
		}
		// Compare data 
		System.out.println("checking all data beans");
		Assert.assertEquals(sourceManager.databeans().size(), manager.databeans().size());
		Assert.assertTrue(compareDataItemsRecursively(sourceManager.getRootFolder(), manager.getRootFolder()));
		
		return true;
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

	
	private boolean compareDataItemsRecursively (DataItem item1, DataItem item2) throws IOException {
		System.out.println("comparing " + item1.getName() + " and " + item2.getName());
		
		// Check name
		Assert.assertEquals(item1.getName(), item2.getName());
		
		// Check that same class
		DataFolder folder1 = null, folder2 = null;
		DataBean bean1 = null, bean2 = null;
		if (item1 instanceof DataFolder) {
			folder1 = (DataFolder) item1;
			Assert.assertTrue(item2 instanceof DataFolder);
			folder2 = (DataFolder) item2;
		} else {
			bean1 = (DataBean) item1;
			Assert.assertTrue(item2 instanceof DataBean);
			bean2 = (DataBean) item2;
		}
	
		// DataBean
		if (item1 instanceof DataBean) {
			// name
			Assert.assertEquals(bean1.getName(),  bean2.getName()); 
			
			// exact size
			if (CHECK_EXACT_OUTPUT_SIZE) {
				Assert.assertEquals(bean1.getContentLength(), bean2.getContentLength());
			}
			
			// content
			if (bean1.getContentLength() > 0) {
				Assert.assertTrue(bean2.getContentLength() > 0);
			}
			
			
//			InputStream in1 = null;
//			InputStream in2 = null;
//
//			try {
//				in1 = bean1.getContentByteStream();
//				in2 = bean2.getContentByteStream();
//				Assert.assertTrue(IOUtils.contentEquals(in1, in2));
//			} finally {
//				IOUtils.closeIfPossible(in1);
//				IOUtils.closeIfPossible(in2);
//			}
//			
			
			
		}
		
		// DataFolder 
		else {
			
			// child count
			
			Assert.assertEquals(folder1.getChildCount(), folder2.getChildCount());
			
			// children equal
			// TODO sort children first?
			// TODO check contents, links
			Iterator<DataItem> iterator1 = folder1.getChildren().iterator();
			Iterator<DataItem> iterator2 = folder2.getChildren().iterator();
			while (iterator1.hasNext()) {
				Assert.assertTrue(compareDataItemsRecursively(iterator1.next(), iterator2.next()));
			}
		}
	
		return true;
	
	}
	
	public static void main(String[] args) throws Exception {
		
		DirectoryLayout.initialiseClientLayout("http://chipster.csc.fi/rc/chipster-config.xml");
		
		// initialize
		SessionReplayTest test = null;
		try {
			test = new SessionReplayTest("nagios", "99b44a74-d735-4588-93cd-314fbc203e0a", "http://chipster.csc.fi/rc/chipster-config.xml");
		} catch (Exception e) {
			System.out.println("CHIPSTER TEST: Initializing the test failed");
			System.exit(3);
		}

		// run
		boolean success = false;
		try {
			test.setUp();
			success = test.test();
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("TEST ERROR");
			return;
		} finally {
			test.tearDown();
		}

		if (success) {
			System.out.println("TEST OK");
			System.exit(0);
		} else {
			System.out.println("TEST FAILED");
			System.exit(1);
		}
	
		
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
