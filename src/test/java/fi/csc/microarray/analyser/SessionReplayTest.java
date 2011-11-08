package fi.csc.microarray.analyser;

import java.awt.event.MouseEvent;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import javax.swing.Icon;
import javax.xml.bind.Unmarshaller;
import javax.xml.transform.stream.StreamSource;

import org.testng.Assert;

import de.schlichtherle.truezip.zip.ZipFile;
import fi.csc.microarray.TestConstants;
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
import fi.csc.microarray.client.operation.ToolModule;
import fi.csc.microarray.client.session.NonStoppingValidationEventHandler;
import fi.csc.microarray.client.session.UserSession;
import fi.csc.microarray.client.session.schema.DataType;
import fi.csc.microarray.client.session.schema.OperationType;
import fi.csc.microarray.client.session.schema.SessionType;
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

public class SessionReplayTest extends MessagingTestBase {

	private File[] testSessions = new File[] {
			new File("/home/hupponen/ngs-test.zip")
	};
	
	public SessionReplayTest(String username, String password, String configURL) {
		super(username, password, configURL);
	}

	public State test() throws Exception {
		
		for (File testSession : testSessions) {

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
			sourceManager.loadSession(testSession, false);

			
			// Pick import operations FIXME 
			LinkedList<OperationRecord> importOperationRecords = new LinkedList<OperationRecord>();
			for (DataBean dataBean : sourceManager.databeans()) {
				OperationRecord operationRecord = dataBean.getOperationRecord();
				if (OperationDefinition.IMPORT_DEFINITION_ID.equals(operationRecord.getNameID().getID()) && 
						! importOperationRecords.contains(operationRecord)) {
					// copy imported databean
					DataBean dataBeanCopy = manager.createDataBean(dataBean.getName(), testSession, dataBean.getContentUrl().getRef());
					dataBeanCopy.setOperationRecord(operationRecord);
					manager.getRootFolder().addChild(dataBeanCopy);
					
					importOperationRecords.add(operationRecord);
				}
			}

			LinkedList<OperationRecord> operationRecords = new LinkedList<OperationRecord>();
			for (DataBean dataBean : sourceManager.databeans()) {
				OperationRecord operationRecord = dataBean.getOperationRecord();
				if (!operationRecords.contains(operationRecord)) {
				operationRecords.add(operationRecord);
				}
			}
			
			System.out.println("databeans: " + sourceManager.databeans().size());
			System.out.println("operation records: " + operationRecords.size());
			
			
			System.out.println("jep");



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
					DataBean inputBean = manager.getDataBean(inputRecord.getValue().getName());
					if (inputBean != null) {
						inputBeans.add(inputBean);
					}
				}

				// check if there are enough inputs
				if (operationRecord.getInputs().size() != inputBeans.size()) {
					System.out.println("not enough inputs for " + operationRecord.getFullName() + ", skipping");
					continue;
				}

				System.out.println("creating operation for " + operationRecord.getFullName());

				Operation operation = new Operation(getOperationDefinition(operationRecord.getNameID().getID(), toolModules), inputBeans.toArray(new DataBean[] {}));
				
				// Parameters

				// Set up task
				Task task = executor.createTask(operation);
				
				// Execute the task
				CountDownLatch latch = new CountDownLatch(1);
				task.addTaskEventListener(new JobResultListener(latch));
				executor.startExecuting(task);
				latch.await(TestConstants.TIMEOUT_AFTER, TimeUnit.MILLISECONDS);
				State endState = task.getState();
				Assert.assertEquals(endState, State.COMPLETED);
				if (!task.getState().equals(State.COMPLETED)) {
					return task.getState();
				}
				
				break;
				// Compare outputs to source session (or should we do this in a one go at the end???)
				// FIXME
			}
		}
		
		return State.COMPLETED;
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
		State state = null;
		try {
			test.setUp();
			state = test.test();

		} catch (Throwable e) {
			try {
				test.tearDown();
			} catch (Exception e2) {
			}
			
			e.printStackTrace();
			String errorMessage = e.getMessage();
			if (errorMessage != null) {
				errorMessage = errorMessage.replace("<", " ");
				errorMessage = errorMessage.replace(">", " ");
				
			}
			
			System.out.println("CHIPSTER FAILED: " + errorMessage);
			System.exit(2);
		}
		try {
			test.tearDown();
		} catch (Exception e) {
		}
	
		if (state.equals(State.COMPLETED)) {
			System.out.println("CHIPSTER OK");
			System.exit(0);
		} else {
			System.out.println("CHIPSTER FAILED: " + "Task COMPLETED, output file did not match");
			System.exit(2);
		}
	}
	
	public static class SessionLoadingSkeletonApplication extends ClientApplication {

		private SessionReplayTest parent;

		public SessionLoadingSkeletonApplication(SessionReplayTest parent, LinkedList<ToolModule> toolModules) {
			this.parent = parent;
			this.toolModules = toolModules;
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
			throw new UnsupportedOperationException("not supported by skeleton app");
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
