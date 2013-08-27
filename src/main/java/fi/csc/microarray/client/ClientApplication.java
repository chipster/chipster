/*
 * Created on Mar 2, 2005
 *
 */
package fi.csc.microarray.client;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseEvent;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import javax.swing.Icon;
import javax.swing.Timer;

import org.apache.log4j.Logger;
import org.eclipse.jetty.util.IO;

import fi.csc.microarray.client.dataimport.ImportItem;
import fi.csc.microarray.client.dataimport.ImportSession;
import fi.csc.microarray.client.dataimport.ImportUtils;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.ChipsterDialog.PluginButton;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.Operation.DataBinding;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.operation.ToolCategory;
import fi.csc.microarray.client.operation.ToolModule;
import fi.csc.microarray.client.selection.DataSelectionManager;
import fi.csc.microarray.client.session.UserSession;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.tasks.Task.State;
import fi.csc.microarray.client.tasks.TaskEventListener;
import fi.csc.microarray.client.tasks.TaskException;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.client.visualisation.Visualisation.Variable;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.client.workflow.WorkflowManager;
import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataChangeEvent;
import fi.csc.microarray.databeans.DataChangeListener;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.SourceMessageListener;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.messaging.auth.ClientLoginListener;
import fi.csc.microarray.module.Module;
import fi.csc.microarray.module.ModuleManager;
import fi.csc.microarray.util.Files;


/**
 * This is the logical essence of Chipster client application. It does
 * not tell how client should look or react, but what and how it should
 * do.
 *  
 * @author Aleksi Kallio
 *
 */
public abstract class ClientApplication {

	protected static final String ALIVE_SIGNAL_FILENAME = "i_am_alive";

	protected static final int MEMORY_CHECK_INTERVAL = 2*1000;
	protected static final int SESSION_BACKUP_INTERVAL = 5 * 1000;

	// Logger for this class
	protected static Logger logger;

    // 
	// ABSTRACT INTERFACE
	//
	protected abstract void initialiseGUIThreadSafely(File mostRecentDeadTempDirectory) throws MicroarrayException, IOException;
	protected abstract void taskCountChanged(int newTaskCount, boolean attractAttention);	
	public abstract void reportExceptionThreadSafely(Exception e);
	public abstract void reportException(Exception e);
	public abstract void reportTaskError(Task job) throws MicroarrayException;
	public abstract void importGroup(Collection<ImportItem> datas, String folderName);
	public abstract DataFolder initializeFolderForImport(String folderName);
	public abstract void showSourceFor(String operationName) throws TaskException;
	public abstract void showHistoryScreenFor(DataBean data);
    public abstract void showDetailsFor(DataBean data);
    public abstract void showPopupMenuFor(MouseEvent e, DataItem data);
    public abstract void showPopupMenuFor(MouseEvent e, List<DataItem> datas);
    public abstract void showImportToolFor(File file, String destinationFolder, boolean skipActionChooser);	
    public abstract void visualiseWithBestMethod(FrameType target);
    public abstract void reportInitialisationThreadSafely(String report, boolean newline);
    public abstract Icon getIconFor(DataItem data);
	public abstract void viewHelp(String id);
	public abstract void viewHelpFor(OperationDefinition operationDefinition);
	public abstract void showDialog(String title, String message, String details, Severity severity, boolean modal);
	public abstract void showDialog(String title, String message, String details, Severity severity, boolean modal, DetailsVisibility detailsVisibility, PluginButton button);
	public abstract void showDialog(String title, String message, String details, Severity severity, boolean modal, DetailsVisibility detailsVisibility, PluginButton button, boolean feedBackEnabled);
	public abstract void deleteDatas(DataItem... datas);	
	public abstract void createLink(DataBean source, DataBean target, Link type);
	public abstract void removeLink(DataBean source, DataBean target, Link type);
	public abstract File saveWorkflow();
	public abstract File openWorkflow();
	public abstract void loadSession();
	public abstract void loadSessionFrom(URL url);
	public abstract void loadSessionFrom(File file);
	public abstract void restoreSessionFrom(File file);
	public abstract void saveSession();
	public abstract void runWorkflow(URL workflowScript);
	public abstract void runWorkflow(URL workflowScript, AtEndListener atEndListener);
	public abstract void flipTaskListVisibility(boolean closeIfVisible); // TODO should not be here (GUI related)
	public abstract void setMaximisedVisualisationMode(boolean maximisedVisualisationMode);
	public abstract VisualisationFrameManager getVisualisationFrameManager();
	public abstract void runBlockingTask(String taskName, final Runnable runnable);
	public abstract DataManager getDataManager();

	/**
	 * Method is called periodically to maintain state that cannot be maintained 
	 * in realtime. 
	 */
	public abstract void checkFreeMemory();
	
	// 
	// CONCRETE IMPLEMENTATIONS (SOME PARTIAL)
	//
	
	/**
	 * Listens to jobExecutor's state in general.
	 */
	private PropertyChangeListener jobExecutorChangeListener = new PropertyChangeListener() {	
		public void propertyChange(PropertyChangeEvent evt) {
			taskCountChanged((Integer)evt.getNewValue(), true);
			logger.debug("JobExecutor property changed event: " + evt.getPropertyName() + ": " + (Integer)evt.getNewValue());
		}		
	};
	
	protected String metadata;
	protected CountDownLatch definitionsInitialisedLatch = new CountDownLatch(1);
	
	private boolean eventsEnabled = false;
	private PropertyChangeSupport eventSupport = new PropertyChangeSupport(this);
	
	protected String requestedModule;

	/**
	 * Tool modules contain tool categories, that contain the tools. 
	 */
	protected LinkedList<ToolModule> toolModules = new LinkedList<ToolModule>(); 
	protected WorkflowManager workflowManager;
	protected DataManager manager;
    protected DataSelectionManager selectionManager;
    protected ServiceAccessor serviceAccessor;
	protected TaskExecutor taskExecutor;
	protected boolean isStandalone;
	private AuthenticationRequestListener overridingARL;

	protected boolean unsavedChanges = false;
	protected boolean unbackuppedChanges = false;

	protected File aliveSignalFile;
	private LinkedList<File> deadDirectories = new LinkedList<File>();

    protected ClientConstants clientConstants;
    protected Configuration configuration;

	private String initialisationWarnings = "";
	
	private String announcementText = null;

	public ClientApplication() {
		this(false, null);
	}

	public ClientApplication(boolean isStandalone, AuthenticationRequestListener overridingARL) {
		this.configuration = DirectoryLayout.getInstance().getConfiguration();
		this.clientConstants = new ClientConstants();
		this.serviceAccessor = isStandalone ? new LocalServiceAccessor() : new RemoteServiceAccessor();
		this.isStandalone = isStandalone;
		this.overridingARL = overridingARL;
	}
    
	protected void initialiseApplication() throws MicroarrayException, IOException {
		
		//Executed outside EDT, modification of Swing forbidden
		
		// these had to be delayed as they are not available before loading configuration
		logger = Logger.getLogger(ClientApplication.class);

		try {

			// Fetch announcements
			fetchAnnouncements();
			
			// Initialise modules
			final ModuleManager modules = new ModuleManager(requestedModule);
			Session.getSession().setModuleManager(modules);

			// Initialise workflows
			this.workflowManager = new WorkflowManager(this);

			// Initialise data management
			this.manager = new DataManager();
			Session.getSession().setDataManager(manager);
			modules.plugAll(this.manager, Session.getSession());
			this.selectionManager = new DataSelectionManager(this);
			Session.getSession().setClientApplication(this);
		
			// try to initialise JMS connection (or standalone services)
			logger.debug("Initialise JMS connection.");
			reportInitialisationThreadSafely("Connecting to broker at " + configuration.getString("messaging", "broker-host") + "...", true);
			serviceAccessor.initialise(manager, getAuthenticationRequestListener());
			this.taskExecutor = serviceAccessor.getTaskExecutor();
			Session.getSession().setServiceAccessor(serviceAccessor);
			reportInitialisationThreadSafely(" ok", false);

			// Check services
			reportInitialisationThreadSafely("Checking remote services...", true);
			String status = serviceAccessor.checkRemoteServices();
			if (!ServiceAccessor.ALL_SERVICES_OK.equals(status)) {
				throw new Exception(status);
			}
			reportInitialisationThreadSafely(" ok", false);
			
			// Fetch descriptions from compute server
			reportInitialisationThreadSafely("Fetching analysis descriptions...", true);
			initialisationWarnings += serviceAccessor.fetchDescriptions(modules.getPrimaryModule());
			toolModules.addAll(serviceAccessor.getModules());

			// Add local modules also when in remote mode
			if (!isStandalone) {
				ServiceAccessor localServiceAccessor = new LocalServiceAccessor();
				localServiceAccessor.initialise(manager, null);
				localServiceAccessor.fetchDescriptions(modules.getPrimaryModule());
				toolModules.addAll(localServiceAccessor.getModules());
			}

			// Add internal operation definitions
			ToolCategory internalCategory = new ToolCategory("Internal tools");
			internalCategory.addOperation(OperationDefinition.IMPORT_DEFINITION);
			internalCategory.addOperation(OperationDefinition.CREATE_DEFINITION);
			ToolModule internalModule = new ToolModule("internal");
			internalModule.addHiddenToolCategory(internalCategory);
			toolModules.add(internalModule);

			// Update to splash screen that we have loaded tools
			reportInitialisationThreadSafely(" ok", false);

			// start listening to job events
			taskExecutor.addChangeListener(jobExecutorChangeListener);

			// definitions are now initialised
			definitionsInitialisedLatch.countDown();
			
			reportInitialisationThreadSafely("Checking session backups...", true);
			File mostRecentDeadTempDirectory = checkTempDirectories();
			reportInitialisationThreadSafely(" ok", false);

			// we can initialise graphical parts of the system
			initialiseGUIThreadSafely(mostRecentDeadTempDirectory);

			// Remember changes to confirm close only when necessary and to backup when necessary
			manager.addDataChangeListener(new DataChangeListener() {
				public void dataChanged(DataChangeEvent event) {
					unsavedChanges = true;
					unbackuppedChanges = true;
				}
			});

			// Start checking amount of free memory 
			final Timer memoryCheckTimer = new Timer(MEMORY_CHECK_INTERVAL, new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					ClientApplication.this.checkFreeMemory();
				}
			});
			memoryCheckTimer.setCoalesce(true);
			memoryCheckTimer.setRepeats(true);
			memoryCheckTimer.setInitialDelay(0);
			memoryCheckTimer.start();

			// Start checking if background backup is needed
			aliveSignalFile = new File(manager.getRepository(), "i_am_alive");
			aliveSignalFile.createNewFile();
			aliveSignalFile.deleteOnExit();

			Timer timer = new Timer(SESSION_BACKUP_INTERVAL, new ActionListener() {
				@Override
				public void actionPerformed(ActionEvent e) {
					aliveSignalFile.setLastModified(System.currentTimeMillis()); // touch the file
					if (unbackuppedChanges) {

						File sessionFile = UserSession.findBackupFile(getDataManager().getRepository(), true);
						sessionFile.deleteOnExit();

						try {
							getDataManager().saveLightweightSession(sessionFile);

						} catch (Exception e1) {
							logger.warn(e1); // do not care that much about failing session backups
						}
					}
					unbackuppedChanges = false;
				}
			});

			timer.setCoalesce(true);
			timer.setRepeats(true);
			timer.setInitialDelay(SESSION_BACKUP_INTERVAL);
			timer.start();

		} catch (Exception e) {
			e.printStackTrace();
			throw new MicroarrayException(e);
		}
	}

	/**
	 * Add listener for applications state changes.
	 */
	public void addClientEventListener(PropertyChangeListener listener) {
		eventSupport.addPropertyChangeListener(listener);		
	}

	/**
	 * @see #addClientEventListener(PropertyChangeListener)
	 */
    public void removeClientEventListener(PropertyChangeListener listener) {
        eventSupport.removePropertyChangeListener(listener);       
    }
    
    public DataSelectionManager getSelectionManager() {
    	return selectionManager;
    }
    
    public void selectAllItems(){
		List<DataBean> datas = manager.databeans();
		for (DataBean data : datas) {
			
			selectionManager.selectMultiple(data, this);
			
		}
    }

	public void setVisualisationMethod(VisualisationMethod method, List<Variable> variables, List<DataBean> datas, FrameType target ) {
		fireClientEvent(new VisualisationMethodChangedEvent(this, method, variables, datas, target));
	}
	
	public void setVisualisationMethod(VisualisationMethodChangedEvent e){
		fireClientEvent(e);
	}
	
	public void setEventsEnabled(boolean eventsEnabled) {
		this.eventsEnabled = eventsEnabled;
		taskExecutor.setEventsEnabled(eventsEnabled);			
	}
	
	/**
	 * Renames the given dataset with the given name and updates the change
	 * on screen.
	 * 
	 * @param data Dataset to rename.
	 * @param newName The new name. Must contain at least one character.
	 */
	public void renameDataItem(DataItem data, String newName) {
		data.setName(newName);
	}
	

	public void executeOperation(final Operation operation) {

		// check if guest user
		if (Session.getSession().getUsername() != null && Session.getSession().getUsername().equals(configuration.getString("security", "guest-username"))) {
			showDialog("Running tools is disabled for guest users.", "",
					null, Severity.INFO, true, DetailsVisibility.DETAILS_ALWAYS_HIDDEN, null);
			return;
		}
		
		// check job count
		if (taskExecutor.getRunningTaskCount() >= clientConstants.MAX_JOBS) {
			showDialog("Task not started as there are maximum number of tasks already running.", "You can only run " + clientConstants.MAX_JOBS + " tasks at the same time. Please wait for one of the currently running tasks to finish and try again.",
						null, Severity.INFO, false);
			return;
		}
		
		// start executing the task
		Task task = taskExecutor.createTask(operation);
		task.addTaskEventListener(new TaskEventListener() {
			public void onStateChange(Task job, State oldState, State newState) {
				if (newState.isFinished()) {
					try {
						// FIXME there should be no need to pass the operation as it goes within the task
						onFinishedTask(job, operation);
					} catch (Exception e) {
						reportException(e);
					}
				}
			}
		});

		try {
			taskExecutor.startExecuting(task);
		} catch (TaskException te) {
			reportException(te);
		}
	}
	
	/**
	 * When a job finishes, this is called by the JobEventListener that
	 * monitors the execution. This creates a new dataset out of the
	 * results and inserts it to the data set views.
	 * 
	 * @param task The finished task.
	 * @param oper The finished operation, which in fact is the GUI's
	 * 			   abstraction of the concrete executed job. Operation
	 * 			   has a decisively longer life span than its
	 * 			   corresponding job entity.
	 * @throws MicroarrayException 
	 * @throws IOException 
	 */
	public void onFinishedTask(Task task, Operation oper) throws MicroarrayException, IOException {
		
		LinkedList<DataBean> newBeans = new LinkedList<DataBean>();
		try {

			logger.debug("operation finished, state is " + task.getState());
			
			if (task.getState() == State.CANCELLED) {
				// task cancelled, do nothing
				
			} else if (!task.getState().finishedSuccesfully()) {
				// task unsuccessful, report it
				reportTaskError(task);
				
			} else {
				// task completed, create datasets etc.
				newBeans = new LinkedList<DataBean>();

				// read operated datas
				Module primaryModule = Session.getSession().getPrimaryModule();
				LinkedList<DataBean> sources = new LinkedList<DataBean>();
				for (DataBinding binding : oper.getBindings()) {
					// do not create derivation links for metadata datasets
					// also do not create links for sources without parents
					// this happens when creating the input databean for example
					// for import tasks
					// FIXME should such a source be deleted here?
					if (!primaryModule.isMetadata(binding.getData()) && (binding.getData().getParent() != null)) {
						sources.add(binding.getData());

					}
				}

				// decide output folder
				DataFolder folder = null;
				if (oper.getOutputFolder() != null) {
					folder = oper.getOutputFolder();
				} else if (sources.size() > 0) {
					for (DataBean source : sources) {
						if (source.getParent() != null) {
							folder = source.getParent();
						}
					}
				}
				// use root if no better option 
				if (folder == null) {
					folder = manager.getRootFolder();
				}


				// read outputs and create derivational links for non-metadata beans
				DataBean metadataOutput = null;
				OperationRecord operationRecord = new OperationRecord(oper);
				operationRecord.setSourceCode(task.getSourceCode());
				
				for (String outputName : task.outputNames()) {

					DataBean output = task.getOutput(outputName);
					output.setOperationRecord(operationRecord);


					// set sources
					for (DataBean source : sources) {
						output.addLink(Link.DERIVATION, source);
					}

					// initialise cache
					try {
						output.initialiseStreamStartCache();
					} catch (IOException e) {
						throw new MicroarrayException(e);
					}

					// connect data (events are generated and it becomes visible)
					folder.addChild(output);

					// check if this is metadata
					// for now this must be after folder.addChild(), as type tags are added there
					if (primaryModule.isMetadata(output)) {
						metadataOutput = output;				
					}
					
					newBeans.add(output);
				}

				// link metadata output to other outputs
				if (metadataOutput != null) {
					for (DataBean bean : newBeans) {
						if (bean != metadataOutput) {
							metadataOutput.addLink(Link.ANNOTATION, bean);
						}
					}

					primaryModule.postProcessOutputMetadata(oper, metadataOutput);				
				}

			}			
	
		} finally {
			
			// notify result listener
			if (oper.getResultListener() != null) {
				if (task.getState().finishedSuccesfully()) {
					oper.getResultListener().resultData(newBeans);
				} else {
					oper.getResultListener().noResults();
				}
			}
		}
	}
	
	protected void quit() {		
		logger.debug("quitting client");
		
		try {
			serviceAccessor.close();
		} catch (Exception e) {
			// do nothing
		}
	}
	
	public void fireClientEvent(PropertyChangeEvent event) {
		logger.debug("dispatching event: " + event);
		if (eventsEnabled) {
			eventSupport.firePropertyChange(event);
		}
	}
	
	public List<File> getWorkflows() {
		return workflowManager.getWorkflows();
	}

	public interface SourceCodeListener {
		public void updateSourceCodeAt(int index, String sourceCode);
	}
	
	public void fetchSourceFor(String[] operationIDs, final SourceCodeListener listener) throws MicroarrayException {
		int i = -1;		
		for (String id : operationIDs) {
			i++;
			logger.debug("describe operation " + id);
			if (id == null) {
				listener.updateSourceCodeAt(i, null);
				continue;
			}
			
			SourceMessageListener sourceListener = null;
			try {
				sourceListener = serviceAccessor.retrieveSourceCode(id);
				String source = sourceListener.waitForResponse(60, TimeUnit.SECONDS);
				listener.updateSourceCodeAt(i, source); // source can be null
				
			} catch (Exception e) {
				throw new MicroarrayException(e);
				
			} finally {
				if (sourceListener != null) {
					sourceListener.cleanUp();
				}
			}
			
		}
	}
		
	public void importWholeDirectory(File root) {
		List<File> onlyFiles = new LinkedList<File>();
		
		for (File file : root.listFiles()) {				
			if (file.isFile()) { //not a folder
				onlyFiles.add(file);
			}
		}
		
		ImportSession importSession = new ImportSession(ImportSession.Source.CLIPBOARD, onlyFiles, root.getName(), true);
		ImportUtils.executeImport(importSession);
	}

	/**
	 * 
	 * @param toolId
	 * @return null if operation definition is not found
	 */
	public OperationDefinition getOperationDefinition(String toolId) {
		for (ToolModule module : toolModules) {
			OperationDefinition tool = module.getOperationDefinition(toolId);
			if (tool != null) {
				return tool;
			}
		}
		return null;
	}

	/**
	 * Get OperationDefinition which best matches the given module and category names.
	 * 
	 * Module is matched before category.
	 * 
	 * @param toolId
	 * @param moduleName
	 * @param categoryName
	 * @return null if not found
	 */
	public OperationDefinition getOperationDefinitionBestMatch(String toolId, String moduleName, String categoryName) {
		
		// module match
		ToolModule preferredModule = getModule(moduleName);
		if (preferredModule != null) {
			OperationDefinition preferredTool = preferredModule.getOperationDefinition(toolId, categoryName);
			
			// module and category match
			if (preferredTool != null) {
				return preferredTool;
			} 
			
			// module match, category mismatch
			else {
				preferredTool = preferredModule.getOperationDefinition(toolId);
				if (preferredTool != null) {
					return preferredTool;
				} 
			}
		} 
		
		// module mismatch
		else {
			OperationDefinition toolWithCategoryMismatch = null;
			for (ToolModule module : toolModules) {
				// try to find tool with matching category, return if found
				OperationDefinition tool = module.getOperationDefinition(toolId, categoryName);
				if (tool != null) {
					return tool;
				}

				// try to find tool with mismatching category
				tool = module.getOperationDefinition(toolId);
				if (tool != null) {
					toolWithCategoryMismatch = tool;
				}
			}

			// matching category not found, return with mismatch, may be null
			return toolWithCategoryMismatch;
		}
		return null;
	}

	
	/**
	 * FIXME Better handling for existing file
	 * 
	 * @param data
	 * @param selectedFile
	 */
	protected void exportToFile(final DataBean data, final File selectedFile) {
		runBlockingTask("exporting file", new Runnable() {

			public void run() {
				try {
					File newFile = selectedFile;
					int i = 1;
					while (newFile.exists()) {
						i++;
						String[] parts = Files.parseFilename(selectedFile); 
						newFile = new File(parts[0] + File.separator + parts[1] + "_" + i + "." + parts[2]);
					}

					newFile.createNewFile();		
					FileOutputStream out = new FileOutputStream(newFile);
					IO.copy(data.getContentByteStream(), out);
					out.close();
				} catch (Exception e) {
					throw new RuntimeException(e);
				}
			}
			
		});
		
	}

	public TaskExecutor getTaskExecutor() {
		return this.taskExecutor;
	}
	
	protected AuthenticationRequestListener getAuthenticationRequestListener() {

		AuthenticationRequestListener authenticator;

		if (overridingARL != null) {
			authenticator = overridingARL;
		} else {
			authenticator = new Authenticator();
		}

		authenticator.setLoginListener(new ClientLoginListener() {
			public void firstLogin() {
			}

			public void loginCancelled() {
				System.exit(1);
			}
		});

		return authenticator;
	}
	
	public boolean isStandalone() {
		return this.isStandalone;
	}

	/**
	 * Collects all dead temp directories and returns the most recent
	 * that has a restorable session .
	 */
	protected File checkTempDirectories() throws IOException {

		Iterable<File> tmpDirectories = getDataManager().listAllRepositories();
		File mostRecentDeadSignalFile = null;
		
		for (File directory : tmpDirectories) {

			// Skip current temp directory
			if (directory.equals(getDataManager().getRepository())) {
				continue;
			}			
			
			// Check is it alive, wait until alive file should have been updated
			File aliveSignalFile = new File(directory, ALIVE_SIGNAL_FILENAME);
			long originalLastModified = aliveSignalFile.lastModified();
			boolean unsuitable = false;
			while ((System.currentTimeMillis() - aliveSignalFile.lastModified()) < 2*SESSION_BACKUP_INTERVAL) {			
				
				// Updated less than twice the interval time ago ("not too long ago"), so keep on checking
				// until we see new update that confirms it is alive, or have waited long
				// enough that the time since last update grows larger than twice the interval.
				
				// Check if restorable
				if (UserSession.findBackupFile(directory, false) == null) {
					// Does not have backup file, so not interesting for backup.
					// Should be removed anyway, but removing empty directories is not
					// important enough to warrant the extra waiting that follows next.
					// So we will skip this and if it was dead, it will be anyway 
					// cleaned away in the next client startup.
					
					unsuitable = true;
					break;
				}
				
				// Check if updated
				if (aliveSignalFile.lastModified() != originalLastModified) {
					unsuitable = true;
					break; // we saw an update, it is alive
				}

				// Wait for it to update
				try {
					Thread.sleep(1000); // 1 second
				} catch (InterruptedException e) {
					// ignore
				}
			}

			if (!unsuitable) {
				// It is dead, might be the one that should be recovered, check that
				deadDirectories.add(directory);
				File deadSignalFile = new File(directory, ALIVE_SIGNAL_FILENAME);
				if (UserSession.findBackupFile(directory, false) != null 
						&& (mostRecentDeadSignalFile == null 
						|| mostRecentDeadSignalFile.lastModified() < deadSignalFile.lastModified())) {

					mostRecentDeadSignalFile = deadSignalFile;

				}
			}
		}
		
		return mostRecentDeadSignalFile != null ? mostRecentDeadSignalFile.getParentFile() : null;
	}
	
	public void clearDeadTempDirectories() {
		
		// Try to clear dead temp directories
		try {
			for (File dir : deadDirectories) {
				Files.delTree(dir);
			}
		} catch (Exception e) {
			reportException(e);
		}

		// Remove them from bookkeeping in any case
		deadDirectories.clear();
	}

	private ToolModule getModule(String moduleName) {
		for (ToolModule toolModule : toolModules) {
			if (toolModule.getModuleName().equals(moduleName)) {
				return toolModule;
			}
		}
		return null;
	}
	
	public String getInitialisationWarnings() {
		return initialisationWarnings;
	}

	private void fetchAnnouncements() {
		new Thread(new Runnable() {

			@Override
			public void run() {
				InputStream input = null;
				try {
					
					URL url = new URL("http://chipster.csc.fi/announcements/client.txt");

					URLConnection connection = url.openConnection();
					connection.setUseCaches(false);
					connection.setConnectTimeout(10*1000);
					connection.setReadTimeout(10*1000);
					input = connection.getInputStream();
					announcementText = org.apache.commons.io.IOUtils.toString(input);
				} catch (Exception e) {
					// could fail for many reasons, not critical
				} finally {
					org.apache.commons.io.IOUtils.closeQuietly(input);
				}
			}
		}).start();
	}
	
	public String getAnnouncementText() {
		return this.announcementText;
	}
	
}
