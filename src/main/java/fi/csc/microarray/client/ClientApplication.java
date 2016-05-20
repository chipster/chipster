/*
 * Created on Mar 2, 2005
 *
 */
package fi.csc.microarray.client;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.beans.PropertyChangeSupport;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;
import java.security.NoSuchAlgorithmException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;

import javax.net.ssl.SSLContext;
import javax.net.ssl.SSLParameters;
import javax.swing.Icon;

import org.apache.log4j.Logger;
import org.eclipse.jetty.util.IO;

import fi.csc.microarray.client.dataimport.ImportItem;
import fi.csc.microarray.client.dataimport.ImportSession;
import fi.csc.microarray.client.dataimport.ImportUtils;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.ChipsterDialog.PluginButton;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.Operation.ResultListener;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.operation.ToolCategory;
import fi.csc.microarray.client.operation.ToolModule;
import fi.csc.microarray.client.selection.DataSelectionManager;
import fi.csc.microarray.client.session.SessionManager;
import fi.csc.microarray.client.session.SessionManager.SessionChangedEvent;
import fi.csc.microarray.client.session.SessionManager.SessionManagerCallback;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.tasks.Task.State;
import fi.csc.microarray.client.tasks.TaskEventListener;
import fi.csc.microarray.client.tasks.TaskException;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.client.visualisation.Visualisation.Variable;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.client.workflow.WorkflowManager;
import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.VisualConstants;
import fi.csc.microarray.databeans.ContentType;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.DataNotAvailableHandling;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.DataBean.Traversal;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.HistoryText;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.ChecksumException;
import fi.csc.microarray.filebroker.ChecksumInputStream;
import fi.csc.microarray.messaging.AuthCancelledException;
import fi.csc.microarray.messaging.SourceMessageListener;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.messaging.auth.ClientLoginListener;
import fi.csc.microarray.module.Module;
import fi.csc.microarray.module.ModuleManager;
import fi.csc.microarray.module.chipster.ChipsterInputTypes;
import fi.csc.microarray.module.chipster.MicroarrayModule;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.SwingTools;
import fi.csc.microarray.util.ThreadUtils;


/**
 * This is the logical essence of Chipster client application. It does
 * not tell how client should look or react, but what and how it should
 * do.
 *  
 * @author Aleksi Kallio
 *
 */
public abstract class ClientApplication {

	protected static final int MEMORY_CHECK_INTERVAL = 2*1000;

	// Logger for this class
	protected static Logger logger;
	
    // 
	// ABSTRACT INTERFACE
	//
	public abstract void initialiseGUIThreadSafely(File mostRecentDeadTempDirectory) throws MicroarrayException, IOException;
	public abstract void reportInitialisationThreadSafely(String report, boolean newline);
	public abstract void reportExceptionThreadSafely(Exception e);
	public abstract void reportException(Exception e);
	public abstract void reportTaskError(Task job) throws MicroarrayException;		
	public abstract void showDialog(String title, String message, String details, Severity severity, boolean modal);
	public abstract void showDialog(String title, String message, String details, Severity severity, boolean modal, DetailsVisibility detailsVisibility, PluginButton button);
	public abstract void showDialog(String title, String message, String details, Severity severity, boolean modal, DetailsVisibility detailsVisibility, PluginButton button, boolean feedBackEnabled);	
	public abstract void runBlockingTask(String taskName, final Runnable runnable);

	
	/**
	 * Gets default visualisation method for selected databeans. The method is
	 * selected by following steps:
	 * 
	 * <ol>
	 * <li>If no dataset is selected, return
	 * <code>VisualisationMethod.NONE</code> </li>
	 * <li>If only one dataset is selected, return the default method for the
	 * data </li>
	 * </li>
	 * <li>If multiple datasets are selected, check the best method for each
	 * dataset. If the best method is same for all selected datasets and it can
	 * be used with multiple data, the best method is returned. </li>
	 * <li>If the best method is not same for all of the datas, try to find
	 * just some method which is suitable for all datas and can be used with
	 * multiple datasets. </li>
	 * <li>If there were no method to fill the requirements above, return
	 * <code>VisualisationMethod.NONE</code> </li>
	 * 
	 * @return default visualisation method which is suitable for all selected
	 *         datasets
	 */
	public VisualisationMethod getDefaultVisualisationForSelection() {
		logger.debug("getting default visualisation");

		try {
			List<DataBean> beans = getSelectionManager().getSelectedDataBeans();

			if (beans.size() == 1) {
				return Session.getSession().getVisualisations().getDefaultVisualisationFor(beans.get(0));
			} else if (beans.size() > 1)
				for (VisualisationMethod method : Session.getSession().getVisualisations().getOrderedDefaultCandidates()) {
					if (!method.getHeadlessVisualiser().isForMultipleDatas()) {
						continue;
					}
					if (method.isApplicableTo(beans)) {
						return method;
					}
				}

			return null;

		} catch (Exception e) {
			reportException(e);
			return null;
		}
	}
	
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
	protected ExecutorService backgroundExecutor = Executors.newCachedThreadPool();
	private AuthenticationRequestListener overridingARL;

    protected ClientConstants clientConstants;
    protected Configuration configuration;

	private String initialisationWarnings = "";
	
	private String announcementText = null;
		
	private SessionManager sessionManager;

	public ClientApplication() {
		this(null);
	}

	public ClientApplication(AuthenticationRequestListener overridingARL) {
		this.configuration = DirectoryLayout.getInstance().getConfiguration();
		this.clientConstants = new ClientConstants();
		this.serviceAccessor = new RemoteServiceAccessor();
		this.overridingARL = overridingARL;
	}
    
	public void initialiseApplication(boolean fast) throws MicroarrayException, IOException {
		
		//Executed outside EDT, modification of Swing forbidden
		
		// these had to be delayed as they are not available before loading configuration
		logger = Logger.getLogger(ClientApplication.class);

		try {

			// Fetch announcements
			fetchAnnouncements();
						
			if (requestedModule == null) {
				requestedModule = MicroarrayModule.class.getName();
			}
			
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
			Session.getSession().setServiceAccessor(serviceAccessor);
			reportInitialisationThreadSafely("Connecting to broker at " + configuration.getString("messaging", "broker-host") + "...", false);
			
			try {
				serviceAccessor.initialise(manager, getAuthenticationRequestListener());
			} catch (Exception e) {
				serviceAccessor.close();
				throw e;
			}
			reportInitialisationThreadSafely(" ok", true);
			
			// send first login
			reportInitialisationThreadSafely("Logging in...", false);
			try {
				serviceAccessor.login();
			} catch (AuthCancelledException ace) {
				serviceAccessor.close();
				System.exit(0);
			}			
			reportInitialisationThreadSafely(" ok", true);

			
			this.taskExecutor = serviceAccessor.getTaskExecutor();

			this.sessionManager = new SessionManager(manager, taskExecutor, serviceAccessor.getFileBrokerClient(), new ClientSessionManagerCallback(this));
			

			if (!fast) {
				// Check services
				reportInitialisationThreadSafely("Checking remote services...", false);
				String status = serviceAccessor.checkRemoteServices();
				if (!ServiceAccessor.ALL_SERVICES_OK.equals(status)) {
					throw new Exception(status);
				}
				reportInitialisationThreadSafely(" ok", true);
			}
			
			// Fetch descriptions from compute server
			reportInitialisationThreadSafely("Fetching analysis descriptions...", false);
			initialisationWarnings += serviceAccessor.fetchDescriptions(modules.getPrimaryModule());
			toolModules.addAll(serviceAccessor.getModules());

			// Add local modules also when in remote mode
			if (!isStandalone()) {
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
			reportInitialisationThreadSafely(" ok", true);

			// definitions are now initialised
			definitionsInitialisedLatch.countDown();
			
			File mostRecentDeadTempDirectory = null;
			
			if (!fast) {
				reportInitialisationThreadSafely("Checking session backups...", false);
				mostRecentDeadTempDirectory = sessionManager.checkTempDirectories();
				reportInitialisationThreadSafely(" ok", true);
			}
			// we can initialise graphical parts of the system
			initialiseGUIThreadSafely(mostRecentDeadTempDirectory);		
			
			// disable http cache (only after initialization, because it makes 
			// icon loading from jar much slower (about 18 seconds for icons in VisualConstants) 
			IOUtils.disableHttpCache();
			
		} catch (Exception e) {
			e.printStackTrace();
			throw new MicroarrayException("Startup failed\nDebug info:\n" + getConnectionDebugInfo(), e);
		}
	}
	
	private String getConnectionDebugInfo() {
		String msg = "";
	
		msg += "\nSystem properties\n";
		for (Object key : System.getProperties().keySet()) {
			msg += key + ": \t" + System.getProperty(key.toString()) + "\n";
		}
		
		try {
			SSLParameters sslParams = SSLContext.getDefault().getSupportedSSLParameters();
			
			msg += "\nProtocols\n";
			for (String protocol : sslParams.getProtocols()) {
				msg += protocol + "\n";
			}
			
			msg += "\nCipher suites\n";
			for (String cipher : sslParams.getCipherSuites()) {
				msg += cipher + "\n";
			}
		} catch (NoSuchAlgorithmException e) {
			logger.error("failed to get ssl debug info", e);
			msg += "failed to get ssl debug info\n";
		}
		
		return msg;
	}
	
	public class ClientSessionManagerCallback implements SessionManagerCallback {

		private ClientApplication app;

		public ClientSessionManagerCallback(ClientApplication app) {
			this.app = app;
		}

		@Override
		public void showDialog(String title, String message, String details,
				Severity severity, boolean modal,
				DetailsVisibility detailsVisibility) {
			app.showDialog(title, message, details, severity, modal, detailsVisibility, null);
		}

		@Override
		public void reportException(Exception e) {
			app.reportException(e);
		}

		@Override
		public void sessionChanged(SessionChangedEvent e) {
			app.fireClientEventThreadSafely(e);
		}

		@Override
		public List<OperationRecord> getUnfinishedJobs() {
			List<OperationRecord> unfinishedJobs = new ArrayList<>();			
			for (Task task : Session.getSession().getApplication().getTaskExecutor().getTasks(true, true)) {
				unfinishedJobs.add(task.getOperationRecord());
			}
			return unfinishedJobs;
		}

		@Override
		public void continueJobs(List<OperationRecord> unifinishedJobs) {
			// unfinishedJobs is null in old sessions
			if (unifinishedJobs != null) {
				for (OperationRecord job : unifinishedJobs) {
					continueOperation(job);
				}
			}
		}		
	}

	/**
	 * Only root folder supported in this implementation.
	 * 
	 * @param folderName subclasses may use this to group imported datasets
	 * @return always root folder in this implementation
	 */
	public DataFolder initializeFolderForImport(String folderName) {
		return manager.getRootFolder();
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
	

	public Task executeOperation(final Operation operation) {

		// check if guest user
		if (!operation.getDefinition().isLocal() 
				&& Session.getSession().getUsername() != null 
				&& Session.getSession().getUsername().equals(configuration.getString("security", "guest-username"))) {
			
			showDialog("Running tools is disabled for guest users.", "",
					null, Severity.INFO, true, DetailsVisibility.DETAILS_ALWAYS_HIDDEN, null);
			return null;
		}
		
		// check job count
		if (taskExecutor.getRunningTaskCount() >= clientConstants.MAX_JOBS) {
			showDialog("Task not started as there are maximum number of tasks already running.", "You can only run " + clientConstants.MAX_JOBS + " tasks at the same time. Please wait for one of the currently running tasks to finish and try again.",
						null, Severity.INFO, false);
			return null;
		}
		
		// this operation needs to be stored in a session before the application
		// is closed to be able to receive results later
		sessionManager.setUnsavedChanges();
		
		OperationRecord operationRecord = new OperationRecord(operation);
		
		// start executing the task
		Task task = taskExecutor.createNewTask(operationRecord, operation.getDefinition().isLocal());
		
		task.addTaskEventListener(new TaskEventListener() {
			public void onStateChange(Task job, State oldState, State newState) {
				if (newState.isFinished()) {
					try {
						onFinishedTask(job, operation.getResultListener(), newState, false);
					} catch (Exception e) {
						reportException(e);
					}
				}
			}
		});

		try {
			onNewTask(task, operation);
			
			taskExecutor.startExecuting(task);
		} catch (TaskException | MicroarrayException | IOException te) {
			reportException(te);
		}
		
		return task;
	}
	
	public Task continueOperation(OperationRecord operationRecord) {

		// assume that all continued operations are remote 
		boolean local = false;
		Task task = taskExecutor.createContinuedTask(operationRecord, local);
		
		task.addTaskEventListener(new TaskEventListener() {
			public void onStateChange(Task job, State oldState, State newState) {
				if (newState.isFinished()) {
					try {						
						// result listener is always null for continued tasks
						onFinishedTask(job, null, newState, false);
					} catch (Exception e) {
						reportException(e);
					}
				}
			}
		});

		try {			
			
			taskExecutor.continueExecuting(task);
			
		} catch (TaskException te) {
			reportException(te);
		}
		return task;
	}
	
	public void onNewTask(Task task, Operation oper) throws MicroarrayException, IOException {
		
		Module primaryModule = Session.getSession().getPrimaryModule();
		
		for (DataBean input : task.getInputDataBeans()) {
			if (primaryModule.isMetadata(input)) {				
				primaryModule.preProcessInputMetadata(oper, input);				
			}
		}
	}
	
	/**
	 * When a job finishes, this is called by the JobEventListener that
	 * monitors the execution. This creates a new dataset out of the
	 * results and inserts it to the data set views.
	 * 
	 * @param task The finished task.
	 * @param resultListener will be notified when the task completes. It is 
	 * ignored if the client is restarted before the task is completed.
	 * @param oper The finished operation, which in fact is the GUI's
	 * 			   abstraction of the concrete executed job. Operation
	 * 			   has a decisively longer life span than its
	 * 			   corresponding job entity.
	 * @param newState 
	 * @throws MicroarrayException 
	 * @throws IOException 
	 */
	public void onFinishedTask(final Task task, final ResultListener resultListener, State state, boolean wait) throws MicroarrayException, IOException {	
		
		logger.debug("operation finished, state is " + state);

		if (state == State.CANCELLED) {
			// task cancelled, do nothing
			// notify result listener
			if (resultListener != null) {
				resultListener.noResults();
			}

		} else if (!state.finishedSuccesfully()) {
			// task unsuccessful, report it
			reportTaskError(task);
			// notify result listener
			if (resultListener != null) {
				resultListener.noResults();
			}

		} else {
			// task completed, create datasets etc.

			// read operated datas
			Module primaryModule = Session.getSession().getPrimaryModule();
			final LinkedList<DataBean> sources = new LinkedList<DataBean>();
			for (DataBean bean : task.getInputDataBeans()) {
				// do not create derivation links for metadata datasets
				// also do not create links for sources without parents
				// this happens when creating the input databean for example
				// for import tasks
				// FIXME should such a source be deleted here?
				if (!primaryModule.isMetadata(bean) && (bean.getParent() != null)) {
					sources.add(bean);

				}
			}

			// decide output folder
			final DataFolder folder = manager.getRootFolder();
			
			// read outputs and create derivational links for non-metadata beans
			
			for (DataBean output : task.getOutputs()) {
			
				output.setOperationRecord(task.getOperationRecord());

				// set sources
				for (DataBean source : sources) {
					output.addLink(Link.DERIVATION, source);
				}
			}

			// create outputs and notify result listener
			Future<?> future = backgroundExecutor.submit(new Runnable() {
				@Override
				public void run() {
					createOutputs(task, sources, folder, resultListener);
				}
			});
			
			if (wait) {
				try {
					future.get();
				} catch (InterruptedException | ExecutionException e) {
					throw new MicroarrayException("error when creating outputs", e);
				}
			}
		}			
	}
	
	private void createOutputs(final Task task, LinkedList<DataBean> sources, DataFolder folder, final ResultListener resultListener) {

		// connect data (events are generated and it becomes visible)
		manager.connectChildren(task.getOutputs(), folder);

		ThreadUtils.runInEDT(new Runnable() {
			@Override
			public void run() {
				Module primaryModule = Session.getSession().getPrimaryModule();
				LinkedList<DataBean> newBeans = new LinkedList<DataBean>();
				DataBean metadataOutput = null;
				for (DataBean output : task.getOutputs()) {
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

					try {
						primaryModule.postProcessOutputMetadata(task.getOperationRecord(), metadataOutput);
					} catch (MicroarrayException | IOException e) {
						reportException(e);
					}				
				}

				// notify result listener
				if (resultListener != null) {
					if (task.getState().finishedSuccesfully()) {
						resultListener.resultData(newBeans);
					}
				}			
			}
		});
	}
	public void quit() {		
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
				String source;
				try {
					source = sourceListener.waitForResponse(60, TimeUnit.SECONDS);	
				} catch (AuthCancelledException ace) {
					return;
				}
				
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
		List<Object> onlyFiles = new LinkedList<Object>();
		
		for (File file : root.listFiles()) {				
			if (file.isFile()) { //not a folder
				onlyFiles.add(file);
			}
		}
		
		ImportSession importSession = new ImportSession(ImportSession.Source.FILE, onlyFiles, root.getName(), true);
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
	
	public void exportToFileAndWait(final DataBean data,
			final File selectedFile) {
		try {
			File newFile = selectedFile;
			int i = 1;
			while (newFile.exists()) {
				i++;
				String[] parts = Files.parseFilename(selectedFile); 
				newFile = new File(parts[0] + File.separator + parts[1] + "_" + i + "." + parts[2]);
			}

			newFile.createNewFile();		
			try (ChecksumInputStream in = Session.getSession().getDataManager().getContentStream(data, DataNotAvailableHandling.EXCEPTION_ON_NA)) {
				try (FileOutputStream out = new FileOutputStream(newFile)) {
					IO.copy(in, out);
				}
				manager.setOrVerifyChecksum(data, in.verifyChecksums());
			}
		} catch (ChecksumException e) {
			reportExceptionThreadSafely(new ChecksumException("checksum validation of the exported file failed", e));
		} catch (Exception e) {
			reportExceptionThreadSafely(e);
		}
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
	
	/**
	 * @return true if client is running in standalone mode (no connection to server).
	 */
	public boolean isStandalone() {
		return false; // standalone not supported
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
	
	public DataManager getDataManager() {
		return manager;
	}
	public LinkedList<ToolModule> getToolModules() {
		return toolModules;
	}
	
	public void importGroupAndWait(final Collection<ImportItem> datas,
			final String folderName) {
		DataBean lastGroupMember = null;

		try {

			for (ImportItem item : datas) {

				String dataSetName = item.getInputFilename();
				ContentType contentType = item.getType();
				Object dataSource = item.getInput();


				// Selects folder where data is imported to, or creates a
				// new one
				DataFolder folder = initializeFolderForImport(folderName);

				// create the DataBean
				DataBean data;
				if (dataSource instanceof File) {
					data = manager.createDataBean(dataSetName, (File) dataSource);
					
				} else if (dataSource instanceof URL) {
					data = manager.createDataBean(dataSetName, ((URL) dataSource));
					
				} else {
					throw new RuntimeException("unknown data source type: " + dataSource.getClass().getSimpleName());
				}

				// set the content type
				data.setContentType(contentType);

				// add the operation (all databeans have their own import
				// operation
				// instance, it would be nice if they would be grouped)
				Operation importOperation = new Operation(OperationDefinition.IMPORT_DEFINITION, new DataBean[] { data });
				data.setOperationRecord(new OperationRecord(importOperation));

				// data is ready now, make it visible
				manager.connectChild(data, folder);

				// Create group links only if both datas are raw type
				if (lastGroupMember != null && ChipsterInputTypes.hasRawType(lastGroupMember) && ChipsterInputTypes.hasRawType(data)) {

					DataBean targetData = data;

					// Link new data to all group linked datas of given cell
					for (DataBean sourceData : lastGroupMember.traverseLinks(new Link[] { Link.GROUPING }, Traversal.BIDIRECTIONAL)) {
						logger.debug("Created GROUPING link between " + sourceData.getName() + " and " + targetData.getName());
						createLink(sourceData, targetData, DataBean.Link.GROUPING);
					}

					// Create link to the given cell after looping to avoid
					// link duplication
					createLink(lastGroupMember, targetData, DataBean.Link.GROUPING);
				}

				lastGroupMember = data;

			}

		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}
	public String getHelpUrl(String page) {
		if (!page.startsWith(HelpMapping.MANUAL_ROOT)) {
			return HelpMapping.MANUAL_ROOT + page;
		}
		return page;
	}
	
	public String getHelpFor(OperationDefinition definition) {
		String url = definition.getHelpURL();
	    if (url != null && !url.isEmpty()) {
	        // Link is stored in operation definition
	        url = definition.getHelpURL();
	    } else {
	        // Mostly for microarray
	        // TODO: consider refactoring so that url is stored in definition
	        // and this "else" branch is not needed
	        url = HelpMapping.mapToHelppage(definition);
	    }
	    return url;
	}
	
	public String getHistoryText(DataBean data, boolean title, boolean name, boolean date, boolean versions, boolean oper, boolean code, boolean notes, boolean param) {
		return new HistoryText(data).getHistoryText(title, name, date, versions, oper, code, notes, param);
	}
	
	public Icon getIconFor(DataItem element) {
		if (element instanceof DataFolder) {
			return VisualConstants.getIcon(VisualConstants.ICON_TYPE_FOLDER);
		} else {
			return Session.getSession().getPrimaryModule().getIconFor((DataBean) element);
		}
	}
	
	public void deleteDatasWithoutConfirming(DataItem... datas) {
		
		// check that we have something to delete
		if (datas.length == 0) {
			return; // no selection, do nothing
		}		
		
		// remove all selections
		getSelectionManager().clearAll(true, this);

		// do actual delete
		for (DataItem data : datas) {
			manager.delete(data);
		}
	}

	public void createLink(DataBean source, DataBean target, Link type) {
		source.addLink(type, target);
	}

	public void removeLink(DataBean source, DataBean target, Link type) {
		source.removeLink(type, target);
	}
	
	public List<File> getWorkflows() {
		return workflowManager.getWorkflows();
	}
	
	public void saveWorkflow(File file) throws IOException {
		workflowManager.saveSelectedWorkflow(file);
	}
	
	public void runWorkflowAndWait(URL workflowScript) throws InterruptedException {
		final CountDownLatch latch = new CountDownLatch(1);
		workflowManager.runScript(workflowScript, new AtEndListener() {			
			@Override
			public void atEnd(boolean success) {
				latch.countDown();
			}
		});
		latch.await();
	}
	
	public void runWorkflow(URL workflowScript, boolean runForEach) {
		if (!runForEach) {
			// Run once
			workflowManager.runScript(workflowScript, null);

		} else {
			// Run for every selected data separately
			
			// Store current selection
			List<DataBean> datas = getSelectionManager().getSelectedDataBeans();

			// Select one by one and run workflow
			for (DataBean data : datas) {

				// Need synchronized latch to wait for each workflow execution
				final CountDownLatch latch = new CountDownLatch(1);
				AtEndListener atEndListener = new AtEndListener() {
					@Override
					public void atEnd(boolean success) {
						logger.debug("workflow run for each: at end");
						latch.countDown();
					}
				};

				// Run it
				getSelectionManager().selectSingle(data, this);
				logger.debug("workflow run for each: selected " + getSelectionManager().getSelectedDataBeans().size());
				workflowManager.runScript(workflowScript, atEndListener);
				try {
					latch.await();
				} catch (InterruptedException e) {
					// Ignore
				}
			}

			// Restore original selection
			logger.debug("workflow run for each: restore original selection");
			Collection<DataItem> items = new LinkedList<DataItem>();
			items.addAll(datas);
			getSelectionManager().selectMultiple(items, this);
		}
	}
	
	public SessionManager getSessionManager() {
		return sessionManager;
	}
	
	public void fireClientEventThreadSafely(final PropertyChangeEvent event) {
		SwingTools.runInEventDispatchThread(new Runnable() {			
			@Override
			public void run() {
				fireClientEvent(event);
			}
		});		
	}
}
