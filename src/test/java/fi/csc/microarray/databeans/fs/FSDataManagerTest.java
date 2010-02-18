package fi.csc.microarray.databeans.fs;

import java.awt.event.MouseEvent;
import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.Collection;
import java.util.List;

import javax.jms.JMSException;
import javax.swing.Icon;

import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import fi.csc.microarray.client.AtEndListener;
import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.dataimport.ImportItem;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.operation.OperationCategory;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.tasks.TaskException;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.util.Files;

public class FSDataManagerTest {

	@BeforeSuite
	public void init() throws IOException, IllegalConfigurationException {
		DirectoryLayout.initialiseClientLayout().getConfiguration();
	}
	
	@Test
	public void testFSDataManagerInitialisation() throws IOException {
		new FSDataManager();
	}

	@Test(groups = {"smoke"} )
	public void testSnapshot() throws IOException, MicroarrayException {

		// create original
		DataManager manager1 = new FSDataManager();
		String beanName1 = "My bean.txt";
		String beanName2 = "My other bean.txt";
		String beanName3 = "My other bean.txt"; // test beans with a same name
		DataBean bean1 = manager1.createDataBean(beanName1);
		DataBean bean2 = manager1.createDataBean(beanName2);
		DataBean bean3 = manager1.createDataBean(beanName3);
		bean1.addLink(Link.DERIVATION, bean2);
		manager1.getRootFolder().addChild(bean1);
		manager1.getRootFolder().addChild(bean2);
		manager1.getRootFolder().addChild(bean3);
		
		// save
		File snap = new File("temp-snapshot");
		if (snap.exists()) {
			Files.delTree(snap);
		}		
		int dataCount = manager1.saveSnapshot(snap, new DummyClientApplication());
		
		// check
		Assert.assertTrue(dataCount == 3);
		Assert.assertTrue(snap.exists());
		
		// load
		DataManager manager2 = new FSDataManager();
		manager2.loadSnapshot(snap, manager2.getRootFolder(), new DummyClientApplication());
		
		// check
		DataFolder root1 = manager1.getRootFolder();
		DataFolder root2 = manager2.getRootFolder();
		Assert.assertEquals(root2.toStringRecursively(0), root1.toStringRecursively(0));
		DataBean newBean1 = null;
		for (DataItem item : root2.getChildren()) {
			if (item.getName().equals(beanName1)) {
				newBean1 = (DataBean)item;
				break;
			}
		}
		Assert.assertNotNull(newBean1.getContentByteStream());
		Assert.assertEquals(newBean1.getLinkTargets(Link.DERIVATION).size(), 1);
	}
	
	@Test(groups = {"smoke"} )
	public void testDataBeanCreation() throws IOException, MicroarrayException {
		FSDataManager manager = new FSDataManager();
		FSDataBean bean = manager.createDataBean("samename.txt", new FileInputStream("examples/affy_example.cel"));
	
		InputStream originalData = new BufferedInputStream(new FileInputStream("examples/affy_example.cel"));
		InputStream beanData = new BufferedInputStream(bean.getContentByteStream());
		
		Assert.assertTrue(Files.equalInputStreamContent(originalData, beanData));

		
		// create more databeans to test duplicate names
		manager.createDataBean("samename.txt", new FileInputStream("examples/affy_example.cel"));
		manager.createDataBean("samename.txt", new FileInputStream("examples/affy_example.cel"));
		
	}
	
	private static class DummyClientApplication extends ClientApplication {

		@Override
		public OperationDefinition locateOperationDefinition(String categoryName, String operationName) {
			// dummy implementation
			return new OperationDefinition("name", new OperationCategory("cat. name"), "description", false);
		}
			

		// REST OF THE CLASS IS JUST EMPTY IMPLEMENTATIONS
				
		@Override
		public void createLink(DataBean source, DataBean target, Link type) {
		}

		@Override
		public void flipTaskListVisibility(boolean closeIfVisible) {
		}

		@Override
		protected AuthenticationRequestListener getAuthenticationRequestListener() {
			return null;
		}

		@Override
		public Icon getIconFor(DataItem data) {
			return null;
		}

		@Override
		public void heartBeat() {
		}

		@Override
		public void loadSession() {
		}

		@Override
		public File openWorkflow() {
			return null;
		}

		@Override
		public void removeLink(DataBean source, DataBean target, Link type) {
		}

		@Override
		public void reportException(Exception e) {
		}

		@Override
		public void reportInitialisation(String report, boolean newline) {
		}

		@Override
		public void reportTaskError(Task job) throws MicroarrayException {
		}

		@Override
		public void saveSession() {
		}

		@Override
		public File saveWorkflow() {
			return null;
		}

		@Override
		public void showDetailsFor(DataBean data) {
		}

		@Override
		public void showDialog(String title, String message, String details, Severity severity, boolean modal) {
		}

		@Override
		public void showDialog(String title, String message, String details, Severity severity, boolean modal, DetailsVisibility detailsVisibility) {
		}

		@Override
		public void showHistoryScreenFor(DataBean data) {
		}

		@Override
		public void showImportToolFor(File file, String destinationFolder, boolean skipActionChooser) {
		}

		@Override
		public void showPopupMenuFor(MouseEvent e, DataItem data) {
		}

		@Override
		public void showPopupMenuFor(MouseEvent e, List<DataItem> datas) {
		}

		@Override
		public void showSourceFor(String operationName) throws TaskException {
		}

		@Override
		protected void taskCountChanged(int newTaskCount, boolean attractAttention) {
		}

		@Override
		public void viewHelp(String id) {
		}

		@Override
		public void viewHelpFor(OperationDefinition operationDefinition) {
		}
	

		public void onException(JMSException arg0) {
		}

		public DataManager getDataManager() {
			return null;
		}

		public DataFolder initializeFolderForImport(String string) {
			return null;
		}


		@Override
		public VisualisationFrameManager getVisualisationFrameManager() {
			return null;
		}


		@Override
		public void setMaximisedVisualisationMode(
				boolean maximisedVisualisationMode) {			
		}


		@Override
		public void visualiseWithBestMethod(FrameType target) {
		}


		@Override
		public void importGroup(Collection<ImportItem> datas, String folderName) {
		}


		@Override
		public void runBlockingTask(String taskName, Runnable runnable) {
		}

		@Override
		public void deleteDatas(DataItem... datas) {
		}

		@Override
		public void loadSessionFrom(URL url) {
		}

		@Override
		public void runWorkflow(URL workflowScript) {
		}


		@Override
		public void runWorkflow(URL workflowScript,
				AtEndListener atEndListener) {
			
		}		
	}
}
