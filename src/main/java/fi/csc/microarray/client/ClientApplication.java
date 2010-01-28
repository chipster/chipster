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
import java.net.URL;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.CountDownLatch;

import javax.jms.JMSException;
import javax.swing.Icon;
import javax.swing.Timer;

import org.apache.log4j.Logger;
import org.mortbay.util.IO;

import fi.csc.microarray.analyser.AnalyserServer;
import fi.csc.microarray.client.dataimport.ImportItem;
import fi.csc.microarray.client.dataimport.ImportSession;
import fi.csc.microarray.client.dataimport.ImportUtils;
import fi.csc.microarray.client.dialog.ChipsterDialog.DetailsVisibility;
import fi.csc.microarray.client.dialog.DialogInfo.Severity;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.operation.OperationCategory;
import fi.csc.microarray.client.operation.OperationDefinition;
import fi.csc.microarray.client.operation.OperationGenerator;
import fi.csc.microarray.client.operation.Operation.DataBinding;
import fi.csc.microarray.client.operation.Operation.ResultListener;
import fi.csc.microarray.client.selection.DataSelectionManager;
import fi.csc.microarray.client.tasks.Task;
import fi.csc.microarray.client.tasks.TaskEventListener;
import fi.csc.microarray.client.tasks.TaskException;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.client.tasks.Task.State;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager;
import fi.csc.microarray.client.visualisation.VisualisationMethod;
import fi.csc.microarray.client.visualisation.VisualisationMethodChangedEvent;
import fi.csc.microarray.client.visualisation.Visualisation.Variable;
import fi.csc.microarray.client.visualisation.VisualisationFrameManager.FrameType;
import fi.csc.microarray.client.visualisation.methods.PhenodataEditor;
import fi.csc.microarray.client.wizard.WizardContext;
import fi.csc.microarray.client.workflow.WorkflowManager;
import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.databeans.DataItem;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.DataBean.Link;
import fi.csc.microarray.databeans.features.table.EditableTable;
import fi.csc.microarray.databeans.features.table.TableBeanEditor;
import fi.csc.microarray.databeans.fs.FSDataManager;
import fi.csc.microarray.description.ParsedVVSADL;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.AdminAPI;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.Node;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.module.DefaultModules;
import fi.csc.microarray.module.Modules;
import fi.csc.microarray.module.chipster.ChipsterVVSADLParser;
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.Strings;


/**
 * This is the logical essence of Chipster client application. It does
 * not tell how client should look or react, but what and how it should
 * do.
 *  
 * @author Aleksi Kallio
 *
 */
public abstract class ClientApplication implements Node, WizardContext {
	private static final int HEARTBEAT_DELAY = 2*1000;

	/**
	 * Logger for this class
	 */
	private static Logger logger;

	public static File SNAPSHOT_DIR = null;
	public static File OLD_SNAPSHOT_DIR = null;
	
    // 
	// ABSTRACT INTERFACE
	//
	protected abstract AuthenticationRequestListener getAuthenticationRequestListener();
	public abstract void reportException(Exception e);
	public abstract void reportTaskError(Task job) throws MicroarrayException;
	protected abstract void taskCountChanged(int newTaskCount, boolean attractAttention);	
	public abstract void importGroup(Collection<ImportItem> datas, String folderName);
	public abstract void showSourceFor(String operationName) throws TaskException;
	public abstract void showHistoryScreenFor(DataBean data);
    public abstract void showDetailsFor(DataBean data);
    public abstract void showPopupMenuFor(MouseEvent e, DataItem data);
    public abstract void showPopupMenuFor(MouseEvent e, List<DataItem> datas);
    public abstract void showImportToolFor(File file, String destinationFolder, boolean skipActionChooser);	
    public abstract void visualiseWithBestMethod(FrameType target);
    public abstract void reportInitialisation(String report, boolean newline);
    public abstract Icon getIconFor(DataItem data);
	public abstract void viewHelp(String id);
	public abstract void viewHelpFor(OperationDefinition operationDefinition);
	public abstract void showDialog(String title, String message, String details, Severity severity, boolean modal);
	public abstract void showDialog(String title, String message, String details, Severity severity, boolean modal, DetailsVisibility detailsVisibility);
	public abstract void deleteDatas(DataItem... datas);	
	public abstract void createLink(DataBean source, DataBean target, Link type);
	public abstract void removeLink(DataBean source, DataBean target, Link type);
	public abstract File saveWorkflow();
	public abstract File openWorkflow();
	public abstract void loadSession();
	public abstract void loadSessionFrom(URL url);
	public abstract void saveSession();
	public abstract void runWorkflow(URL workflowScript);
	public abstract void runWorkflow(URL workflowScript, AtEndListener atEndListener);
	public abstract void flipTaskListVisibility(boolean closeIfVisible); // TODO should not be here (GUI related)
	public abstract void setMaximisedVisualisationMode(boolean maximisedVisualisationMode);
	public abstract VisualisationFrameManager getVisualisationFrameManager();
	public abstract void runBlockingTask(String taskName, final Runnable runnable);

	/**
	 * Method is called periodically to maintain state that cannot be maintained 
	 * in realtime. 
	 */
	public abstract void heartBeat();
	
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
	
	private NodeBase nodeSupport = new NodeBase() {
		public String getName() {
			return "client";
		}
	};

	protected Collection<OperationCategory> parsedCategories;
	protected String metadata;
	protected CountDownLatch definitionsInitialisedLatch = new CountDownLatch(1);
	
	private boolean eventsEnabled = false;
	private PropertyChangeSupport eventSupport = new PropertyChangeSupport(this);

	protected WorkflowManager workflowManager;
	protected TaskExecutor taskExecutor;
	protected MessagingEndpoint endpoint;
	protected DataManager manager;
    protected DataSelectionManager selectionManager;

    protected ClientConstants clientConstants;
    protected Configuration configuration;

	public ClientApplication() {
		this.configuration = DirectoryLayout.getInstance().getConfiguration();
		this.clientConstants = new ClientConstants();
	}
    
	protected void initialiseApplication() throws MicroarrayException, IOException {
		
		// these had to be delayed as they are not available before loading configuration
		logger = Logger.getLogger(ClientApplication.class);
		SNAPSHOT_DIR = new File(DirectoryLayout.getInstance().getUserDataDir().getAbsolutePath(), "session-snapshot.zip");
		OLD_SNAPSHOT_DIR = new File(DirectoryLayout.getInstance().getUserDataDir().getAbsolutePath(), "workspace-snapshot");
		
		// initialise modules
		Modules modules = DefaultModules.getDefaultModules();
		Session.getSession().putObject("modules", modules);
		
		// initialise workflows
		this.workflowManager = new WorkflowManager(this);
		 
		// initialise data management
		this.manager = new FSDataManager();
		modules.plugFeatures(this.manager);
		Session.getSession().putObject("data-manager", manager);

        this.selectionManager = new DataSelectionManager(this);
		Session.getSession().putObject("application", this);
		
		try {
			// try to initialise JMS connection
			logger.debug("Initialise JMS connection.");
			reportInitialisation("Connecting to broker at " + configuration.getString("messaging", "broker-host") + "...", true);
			this.endpoint = new MessagingEndpoint(this, getAuthenticationRequestListener());
			reportInitialisation(" connected", false);				
			
			//	put network stuff to session
			Session.getSession().putObject("client-endpoint", endpoint);
			taskExecutor = new TaskExecutor(endpoint, manager);
			Session.getSession().putObject("client-job-executor", taskExecutor);
			
			reportInitialisation("Checking remote services...", true);				
			AdminAPI api = new AdminAPI(endpoint.createTopic(Topics.Name.ADMIN_TOPIC, AccessMode.READ_WRITE), null);
			if (!api.areAllServicesUp(true)) {
				throw new Exception("required services are not available (" + api.getErrorStatus() + ")");
			}				
			reportInitialisation(" all are available", false);
			
			// create metadata fetching job
			reportInitialisation("Fetching analysis descriptions...", true);
			final Task describeOperations = taskExecutor.createTask("describe", true);
			
			// run the job (blocking while it is progressing)
			taskExecutor.execute(describeOperations);
			
			// parse metadata
			DataBean metadataBean = describeOperations.getOutput(AnalyserServer.DESCRIPTION_OUTPUT_NAME);
			this.metadata = new String(metadataBean.getContents());
			manager.delete(metadataBean); // don't leave the bean hanging around
			logger.debug("got metadata: " + this.metadata.substring(0, 50) + "...");
			List<ParsedVVSADL> descriptions = new ChipsterVVSADLParser().parseMultiple(this.metadata);
			this.parsedCategories = new OperationGenerator().generate(descriptions).values();
			
			logger.debug("created " + this.parsedCategories.size() + " operation categories");
			
			reportInitialisation(" received and processed", false);
			definitionsInitialisedLatch.countDown();
			
			// start listening to job events
			taskExecutor.addChangeListener(jobExecutorChangeListener);
			
			// start heartbeat
			final Timer timer = new Timer(HEARTBEAT_DELAY, new ActionListener() {
				public void actionPerformed(ActionEvent e) {
					ClientApplication.this.heartBeat();
				}
			});
			timer.setCoalesce(true);
			timer.setRepeats(true);
			timer.setInitialDelay(0);
			timer.start();
			
		} catch (Exception e) {
			showDialog("Starting Chipster failed.", "There could be a problem with the network connection, or the remote services could be down. " +
					"Please see the details below for more information about the problem.\n\n" + 
					"Chipster also fails to start if there has been a version update with a change in configurations. In such case please delete Chipster application settings directory.",
					e.toString(), Severity.ERROR, false);
			
			throw new MicroarrayException(e);
		}


	}
	
	/**
	 * Add listener for applications state changes.
	 */
	public void addPropertyChangeListener(PropertyChangeListener listener) {
		eventSupport.addPropertyChangeListener(listener);		
	}

	/**
	 * @see #addPropertyChangeListener(PropertyChangeListener)
	 */
    public void removePropertyChangeListener(PropertyChangeListener listener) {
        eventSupport.removePropertyChangeListener(listener);       
    }
    
    public List<DataBean> getAllDataBeans(){
		List<DataBean> datas = new ArrayList<DataBean>();
		// The depth of the file structure is max 2, so we don't need recursion
		
		// Iterate the folders
		for (DataItem item1 : this.manager.getRootFolder().getChildren()) {
			if(item1 instanceof DataFolder){
				DataFolder folder = (DataFolder)item1;
				
				// Iterate the datas
				for(DataItem item2 : folder.getChildren()){
					if(item2 instanceof DataBean){
						DataBean bean = (DataBean)item2;
						datas.add(bean);
					}
				}
			}
		}
		return datas;
    }

    public DataSelectionManager getSelectionManager() {
    	return selectionManager;
    }
    
    public void selectAllItems(){
		List<DataBean> datas = getAllDataBeans();
		for (DataBean data : datas) {
			
			selectionManager.selectMultiple(data, this);
			
		}
    }

	public void setVisualisationMethod(VisualisationMethod method, List<Variable> variables, List<DataBean> datas, FrameType target ) {
		dispatchVisualisationEvent(new VisualisationMethodChangedEvent(this, method, variables, datas, target));
	}
	
	public void setVisualisationMethod(VisualisationMethodChangedEvent e){
		dispatchEvent(e);
	}
	
	public void setEventsEnabled(boolean eventsEnabled) {
		this.eventsEnabled = eventsEnabled;
		taskExecutor.setEventsEnabled(eventsEnabled);			
	}
	
	public String getName() {
		return nodeSupport.getName();
	}
	
	public String getHost() {
		return nodeSupport.getHost();
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
	
	public void executeOperation(final OperationDefinition operationDefinition, ResultListener resultListener) {
		
		try {
			Operation operation = new Operation(operationDefinition, getSelectionManager().getSelectedDatasAsArray());
			operation.setResultListener(resultListener);
			executeOperation(operation);
			
		} catch (MicroarrayException e) {
			reportException(e);
		}
	}
	
	public void executeOperation(final Operation operation) {

		// check operation (relevant only for workflows)
		if (operation.getBindings() == null) {
			throw new RuntimeException("Attempted to run " + operation.getDefinition().getFullName() + " with input datasets that were not compatitible with the operation.");
		}
		
		// check job count
		if (taskExecutor.getRunningTaskCount() >= clientConstants.MAX_JOBS) {
			showDialog("Task not started as there are maximum number of tasks already running.", "You can only run " + clientConstants.MAX_JOBS + " tasks at the same time. Please wait for one of the currently running tasks to finish and try again.",
						null, Severity.INFO, false);
			return;
		}
		
		// check job size
		long bytes = 0;
		for (DataBinding binding : operation.getBindings()) {
			bytes += binding.getData().getContentLength();			
		}
		int megabytes = (int)(bytes/1000000L);
		if (megabytes > clientConstants.MAX_JOB_SIZE_MB) {
			showDialog("Task not started since input datasets are too large.", "Maximum size for input datasets is " + clientConstants.MAX_JOB_SIZE_MB + " megabytes.", "Input datasets size: " + megabytes, Severity.INFO, false);
			return;
		}
		
		// execute the job
		operation.execute(new TaskEventListener() {
			public void onStateChange(Task job, State oldState, State newState) {
				if (newState.isFinished()) {
					try {
						onFinishedTask(job, operation);
					} catch (Exception e) {
						reportException(e);
					}
				}
			}
		});				
	}
	
	/**
	 * When a job finishes, this is called by the JobEventListener that
	 * monitors the execution. This creates a new dataset out of the
	 * results and inserts it to the data set views.
	 * 
	 * @param job The finished job.
	 * @param oper The finished operation, which in fact is the GUI's
	 * 			   abstraction of the concrete executed job. Operation
	 * 			   has a decisively longer life span than its
	 * 			   corresponding job entity.
	 * @throws MicroarrayException 
	 * @throws IOException 
	 */
	public void onFinishedTask(Task job, Operation oper) throws MicroarrayException, IOException {
		
		LinkedList<DataBean> newBeans = new LinkedList<DataBean>();
		try {

			logger.debug("operation finished, state is " + job.getState());
			
			// for canceled tasks, do nothing
			if (job.getState() == State.CANCELLED) {

			}
			// for unsuccessful tasks, report failing
			else if (!job.getState().finishedSuccesfully()) { 
				reportTaskError(job);
			}

			// for completed tasks, create datasets etc.
			else {

				newBeans = new LinkedList<DataBean>();

				// read operated datas
				LinkedList<DataBean> sources = new LinkedList<DataBean>();
				for (DataBinding binding : oper.getBindings()) {
					// remove derivation links that start from phenodata
					if (!binding.getData().queryFeatures("/phenodata").exists()) {
						sources.add(binding.getData());

					}
				}

				// decide output folder
				DataFolder folder;
				if (oper.getOutputFolder() != null) {
					folder = oper.getOutputFolder();
				} else if (sources.size() > 0) {
					folder = sources.get(0).getParent();
				} else {
					folder = manager.getRootFolder();
				}


				DataBean phenodata = null;

				for (String outputName : job.outputNames()) {

					DataBean result = job.getOutput(outputName);
					result.setOperation(oper);

					if (result.queryFeatures("/phenodata").exists()) {
						phenodata = job.getOutput(outputName);					
					}

					// set sources
					for (DataBean source : sources) {
						result.addLink(Link.DERIVATION, source);
					}

					// initialise cache
					try {
						result.initialiseStreamStartCache();
					} catch (IOException e) {
						throw new MicroarrayException(e);
					}

					// connect data (events are generated and it becomes visible)
					folder.addChild(result);

					newBeans.add(result);
				}

				if (phenodata != null) {
					// link phenodata to other datasets
					for (DataBean bean : newBeans) {
						if (bean != phenodata) {
							phenodata.addLink(Link.ANNOTATION, bean);
						}
					}

					// if original names are not already contained in the phenodata 
					if (!phenodata.queryFeatures("/column/" + PhenodataEditor.PHENODATA_NAME_COLUMN).exists()) {
						// augment phenodata with original dataset names (using parameter bindings)
						HashSet<String> insertedNames = new HashSet<String>();
						TableBeanEditor tableEditor = new TableBeanEditor(phenodata);
						EditableTable editableTable = tableEditor.getEditable();
						LinkedList<String> newColumn = new LinkedList<String>();
						newColumn.addAll(Arrays.asList(Strings.repeatToArray("", editableTable.getRowCount())));
						editableTable.addColumn(PhenodataEditor.PHENODATA_NAME_COLUMN, 1, newColumn); // add after sample column 
						for (int ri = 0; ri < editableTable.getRowCount(); ri++) {
							String sample = editableTable.getValue(PhenodataEditor.PHENODATA_SAMPLE_COLUMN, ri);
							boolean correctRowFound = false;
							String originalName = null;
							for (DataBinding binding : oper.getBindings()) {
								if (binding.getName().equals(sample)) {
									originalName = binding.getData().getName();
									correctRowFound = true;
									break;
								}
							}
							if (!correctRowFound) {
								originalName = sample; // just duplicate the sample name if proper is not found
							}

							// check that original names are unique
							if (insertedNames.contains(originalName)) {
								final String separator = "/";
								int i = 2;
								while (insertedNames.contains(originalName + separator + i)) {
									i++;
								}
								originalName = originalName + separator + i;
							}

							editableTable.setValue(PhenodataEditor.PHENODATA_NAME_COLUMN, ri, originalName);
							insertedNames.add(originalName);

						}
						tableEditor.write();
					}

					// if chip descriptions (visualisation view names) aren't there already 
					if (!phenodata.queryFeatures("/column/" + PhenodataEditor.PHENODATA_DESCRIPTION_COLUMN).exists()) {
						// copy original dataset names
						TableBeanEditor tableEditor = new TableBeanEditor(phenodata);
						EditableTable editableTable = tableEditor.getEditable();
						LinkedList<String> newColumn = new LinkedList<String>();
						newColumn.addAll(Arrays.asList(Strings.repeatToArray("", editableTable.getRowCount())));
						editableTable.addColumn(PhenodataEditor.PHENODATA_DESCRIPTION_COLUMN, newColumn); 
						for (int ri = 0; ri < editableTable.getRowCount(); ri++) {
							String sample = editableTable.getValue(PhenodataEditor.PHENODATA_NAME_COLUMN, ri);										
							editableTable.setValue(PhenodataEditor.PHENODATA_DESCRIPTION_COLUMN, ri, sample);
						}
						tableEditor.write();
					}				
				}

			}			
	
		} finally {
			if (oper.getResultListener() != null) {
				if (job.getState().finishedSuccesfully()) {
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
			endpoint.close();
		} catch (JMSException je) {
			// do nothing
		}
	}
	
	public void dispatchVisualisationEvent(VisualisationMethodChangedEvent event) {
		logger.debug("VisualisationEvent dispatched: " + event.getNewMethod());
		this.dispatchEvent(event);
	}
	
	public void dispatchEvent(PropertyChangeEvent event) {
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
	
	public void fetchSourceFor(String[] operationNames, final SourceCodeListener listener) throws MicroarrayException {
		try {
			int i = -1;		
			for (String name : operationNames) {
				i++;
				logger.debug("describe operation " + name);
				if (name == null) {
					listener.updateSourceCodeAt(i, null);
					continue;
				}
				final Task describeTask = taskExecutor.createTask("describe-operation", true);
				final int index = i;
				describeTask.addParameter("name", name);
				describeTask.addTaskEventListener(new TaskEventListener() {
					public void onStateChange(Task job, State oldState, State newState) {
						if (newState == State.COMPLETED) {
							try {
								DataBean sourceBean = describeTask.getOutput(AnalyserServer.SOURCECODE_OUTPUT_NAME);
								String source = new String(sourceBean.getContents());
								manager.delete(sourceBean); // don't leave it hanging around
								logger.debug(source);
								listener.updateSourceCodeAt(index, source);
							} catch (MicroarrayException e) {
								reportException(e);
							}
						}
					}
				});
				taskExecutor.startExecuting(describeTask);
			}
		} catch (TaskException e) {
			throw new MicroarrayException(e);
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

	public OperationDefinition locateOperationDefinition(String categoryName, String operationName) {
		for (OperationCategory category : parsedCategories) {
			if (category.getName().equals(categoryName)) {
				for (OperationDefinition definition : category.getOperationList()) {
					if (definition.getName().equals(operationName)) {
						return definition;
					}
				}
			}
		}
		return null;
	}

	public void loadOldSnapshot() throws IOException, MicroarrayException {
		manager.loadOldSnapshot(OLD_SNAPSHOT_DIR, manager.getRootFolder(), this);
	}

	public Iterable<OperationDefinition> getOperationDefinitions() {
		LinkedList<OperationDefinition> definitions = new LinkedList<OperationDefinition>();
		
		for (OperationCategory category: parsedCategories) {
			for (OperationDefinition operationDefinition: category.getOperationList()) {
				definitions.add(operationDefinition);
			}
		}
		return definitions;
	}
	
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
					throw new RuntimeException();
				}
			}
			
		});
		
	}

	public TaskExecutor getTaskExecutor() {
		return this.taskExecutor;
	}
	
}