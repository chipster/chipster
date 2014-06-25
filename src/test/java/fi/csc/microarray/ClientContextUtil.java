package fi.csc.microarray;

import java.io.File;
import java.io.IOException;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.ChipsterDialog.PluginButton;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.operation.ToolCategory;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.ContentType;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.module.ModuleManager;

public class ClientContextUtil {

	private static class SkeletonApplication extends ClientApplication {

		@Override
		public void initialiseGUIThreadSafely(File backupSession) throws MicroarrayException, IOException {
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
		public void reportInitialisationThreadSafely(String report, boolean newline) {
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
		public void runBlockingTask(String taskName, Runnable runnable) {
			// do nothing
		}

		@Override
		public DataManager getDataManager() {
			return null;
		}

		@Override
		public void reportExceptionThreadSafely(Exception e) {
			// TODO Auto-generated method stub
			
		}
	}
	

	/**
	 * Pushes skeleton implementations of different client side classes to
	 * Session, to create a context where client side code can be tested.
	 */
	public static void setupClientContext() throws Exception {
		DirectoryLayout.uninitialise();
		DirectoryLayout.initialiseSimpleLayout().getConfiguration();			
		DataManager manager = new DataManager();
		Session.getSession().setDataManager(manager);
		ModuleManager moduleManager = new ModuleManager("fi.csc.microarray.module.chipster.MicroarrayModule");
		moduleManager.plugAll(manager, Session.getSession());
		Session.getSession().setModuleManager(moduleManager);
		Session.getSession().setClientApplication(new SkeletonApplication());

	}


	/**
	 * Initialise databean metadata etc. with empty values.
	 */
	public static void setupDatabean(DataBean data) throws MicroarrayException {
		data.setContentType(new ContentType("", false, false, "", null, ""));
		OperationRecord operationRecord = new OperationRecord(new Operation(new OperationDefinition("", "", new ToolCategory(""), "", false), new DataBean[] { data } ));
		data.setOperationRecord(operationRecord);
		operationRecord.setModule("");

		
	}
}
