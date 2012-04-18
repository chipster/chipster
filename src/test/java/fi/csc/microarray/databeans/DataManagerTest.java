package fi.csc.microarray.databeans;

import java.awt.event.MouseEvent;
import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.Collection;
import java.util.List;

import javax.swing.Icon;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import fi.csc.microarray.client.AtEndListener;
import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dataimport.ImportItem;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.ChipsterDialog.PluginButton;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.operation.ToolCategory;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.tasks.TaskException;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.ModuleManager;

public class DataManagerTest {

	private DataManager manager; 
	
	private static class SessionTestingSkeletonApplication extends ClientApplication {

		@Override
		protected void initialiseGUI() throws MicroarrayException, IOException {
			// do nothing
		}

		@Override
		protected void taskCountChanged(int newTaskCount, boolean attractAttention) {
			// do nothing
		}

		@Override
		public void reportException(Exception e) {
			// do nothing
		}

		@Override
		public void reportTaskError(Task job) throws MicroarrayException {
			// do nothing
		}

		@Override
		public void importGroup(Collection<ImportItem> datas, String folderName) {
			// do nothing
		}

		@Override
		public DataFolder initializeFolderForImport(String folderName) {
			return null;
		}

		@Override
		public void showSourceFor(String operationName) throws TaskException {
			// do nothing
		}

		@Override
		public void showHistoryScreenFor(DataBean data) {
			// do nothing
		}

		@Override
		public void showDetailsFor(DataBean data) {
			// do nothing
		}

		@Override
		public void showPopupMenuFor(MouseEvent e, DataItem data) {
			// do nothing
		}

		@Override
		public void showPopupMenuFor(MouseEvent e, List<DataItem> datas) {
			// do nothing
		}

		@Override
		public void showImportToolFor(File file, String destinationFolder, boolean skipActionChooser) {
			// do nothing
		}

		@Override
		public void visualiseWithBestMethod(FrameType target) {
			// do nothing
		}

		@Override
		public void reportInitialisation(String report, boolean newline) {
			// do nothing
		}

		@Override
		public Icon getIconFor(DataItem data) {
			return null;
		}

		@Override
		public void viewHelp(String id) {
			// do nothing
		}

		@Override
		public void viewHelpFor(OperationDefinition operationDefinition) {
			// do nothing
		}

		@Override
		public void showDialog(String title, String message, String details, Severity severity, boolean modal) {
			// do nothing
		}

		@Override
		public void showDialog(String title, String message, String details, Severity severity, boolean modal, DetailsVisibility detailsVisibility, PluginButton button) {
			// do nothing
		}

		@Override
		public void showDialog(String title, String message, String details, Severity severity, boolean modal, DetailsVisibility detailsVisibility, PluginButton button, boolean feedBackEnabled) {
			// do nothing
		}

		@Override
		public void deleteDatas(DataItem... datas) {
			// do nothing
		}

		@Override
		public void createLink(DataBean source, DataBean target, Link type) {
			// do nothing
		}

		@Override
		public void removeLink(DataBean source, DataBean target, Link type) {
			// do nothing
		}

		@Override
		public File saveWorkflow() {
			return null;
		}

		@Override
		public File openWorkflow() {
			return null;
		}

		@Override
		public void loadSession() {
			// do nothing
		}

		@Override
		public void loadSessionFrom(URL url) {
			// do nothing
		}

		@Override
		public void restoreSessionFrom(File file) {
			// do nothing
		}

		@Override
		public void saveSession(boolean quit, boolean lightweight) {
			// do nothing
		}

		@Override
		public void runWorkflow(URL workflowScript) {
			// do nothing
		}

		@Override
		public void runWorkflow(URL workflowScript, AtEndListener atEndListener) {
			// do nothing
		}

		@Override
		public void flipTaskListVisibility(boolean closeIfVisible) {
			// do nothing
		}

		@Override
		public void setMaximisedVisualisationMode(boolean maximisedVisualisationMode) {
			// do nothing
		}

		@Override
		public VisualisationFrameManager getVisualisationFrameManager() {
			return null;
		}

		@Override
		public void runBlockingTask(String taskName, Runnable runnable) {
			// do nothing
		}

		@Override
		public DataManager getDataManager() {
			// TODO Auto-generated method stub
			return null;
		}

		@Override
		public void checkFreeMemory() {
			// do nothing
		}
		
	}
	
	@BeforeClass(alwaysRun = true)
	public void init() throws IOException, IllegalConfigurationException, InstantiationException, IllegalAccessException, ClassNotFoundException {
		DirectoryLayout.initialiseSimpleLayout().getConfiguration();			
		this.manager = new DataManager();
		ModuleManager moduleManager = new ModuleManager("fi.csc.microarray.module.chipster.MicroarrayModule");
		moduleManager.plugAll(this.manager, null);
		Session.getSession().setModuleManager(moduleManager);
		Session.getSession().setClientApplication(new SessionTestingSkeletonApplication());
	}
	
	@Test(groups = {"unit"} )
	public void testDataTypes() throws IOException {
		File file = File.createTempFile("test", ".png");
		Assert.assertEquals(manager.guessContentType(file).getType(), "image/png");
	}

	
	@Test(groups = {"unit"} )
	public void testRemoteSessions() throws Exception {
		
		// populate with crap
		File content = File.createTempFile("content", ".tsv");
		DataBean data = manager.createDataBean("test-content", content);
		OperationRecord operationRecord = new OperationRecord(new Operation(new OperationDefinition("", "", new ToolCategory(""), "", false), new DataBean[] { data } ));
		data.setOperationRecord(operationRecord);
		operationRecord.setModule("");
		manager.getRootFolder().addChild(data);
		
		// save
		File session = File.createTempFile("test-remote-session", ".zip");
		manager.saveLightweightSession(session);

		// clear
		manager.deleteAllDataItems();

		// load
		manager.loadSession(session, true);
		
		// check
		Assert.assertEquals(manager.getRootFolder().getChildCount(), 1);
		Assert.assertEquals(manager.getRootFolder().getChildren().iterator().next().getName(), "test-content");
		
		// clean up
		manager.deleteAllDataItems();
	}

}
