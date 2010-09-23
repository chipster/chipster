package fi.csc.microarray.client.tasks;

import java.io.File;

import javax.jms.JMSException;

import fi.csc.chipster.tools.gbrowser.TsvSorter;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.tasks.Task.State;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ElandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileDefinition;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.messaging.TempTopicMessagingListener;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.util.IOUtils.CopyProgressListener;

public class LocalTaskExecutor extends TaskExecutor {

	public LocalTaskExecutor(DataManager manager) throws JMSException {
		super(manager);
	}

	@Override
	public void kill(Task task) {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public void killAll() {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public void startExecuting(final Task task) throws TaskException {
		if (!task.getOperationID().equals("PreprocessNGSSingle.java")) {
			throw new UnsupportedOperationException();
		}

//		task.setState(State.RUNNING);

		Runnable taskRunnable = new Runnable() {
			public void run() {
				DataManager dataManager = Session.getSession().getDataManager();
				ElandParser parser = new ElandParser();
				FileDefinition fileDefinition = parser.getFileDefinition();
				
				try {

					for (DataBean inputDataBean: task.getInputs()) {
						File inputFile = dataManager.getLocalFile(inputDataBean);
						String outputName; 
						int fileExtensionStartPosition = inputFile.getName().lastIndexOf(".");
						if (fileExtensionStartPosition != -1) {
							outputName = inputFile.getName().substring(0, fileExtensionStartPosition) + "-preprocessed" + 
											inputFile.getName().substring(fileExtensionStartPosition, inputFile.getName().length());
						} else {
							outputName = inputFile.getName() + "-preprocessed";
						}
						File outputFile = dataManager.createNewRepositoryFile(outputName);		

						// run sorter
						new TsvSorter().sort(inputFile, outputFile,	fileDefinition.indexOf(ColumnType.CHROMOSOME), fileDefinition.indexOf(ColumnType.BP_START));

						DataBean outputBean = dataManager.createDataBean(outputName, outputFile);
						outputBean.setOperation(new Operation(Session.getSession().getApplication().getOperationDefinition("PreprocessNGSSingle.java"), new DataBean[] {}));
						dataManager.getRootFolder().addChild(outputBean);
						
					}
				} catch (Exception e) {
					throw new RuntimeException(e);
				}

//				task.setState(State.COMPLETED);

			}
		};
	
		
		Session.getSession().getApplication().runBlockingTask("running " + task.getNamePrettyPrinted(), taskRunnable);
		
		
	}
	
	@Override
	public void startExecuting(Task task, int timeout) throws TaskException {
		throw new UnsupportedOperationException();
	}
}
