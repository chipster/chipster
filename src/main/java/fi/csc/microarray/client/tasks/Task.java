/*
 * Created on Feb 10, 2005
 *
 */
package fi.csc.microarray.client.tasks;

import java.util.ArrayList;
import java.util.LinkedList;
import java.util.List;
import java.util.UUID;

import javax.swing.SwingUtilities;

import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.operation.OperationRecord.ParameterRecord;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;

/**
 * @author Aleksi Kallio, Taavi Hupponen
 *
 */
public class Task {

	public enum State {
		NEW("New"), 
		TRANSFERRING_INPUTS("Transferring inputs"),
		WAITING("Waiting"),
		RUNNING("Running"), 
		TRANSFERRING_OUTPUTS("Transferring outputs"),
		COMPLETED("Completed"), 
		FAILED("Failed"),
		FAILED_USER_ERROR("Failed, fixable"),
		ERROR("Error"),
		CANCELLED("Cancelled"), 
		TIMEOUT("Timeout");
		
		private State(String description) {
			this.description = description;
		}
		
		private String description;
		
		public boolean isFinished() {
			if (this.equals(COMPLETED) || this.equals(FAILED) || this.equals(FAILED_USER_ERROR) ||
					this.equals(ERROR) || this.equals(CANCELLED) || this.equals(TIMEOUT)) {
				return true;
			} else {
				return false;
			}
		}
		
		public boolean finishedSuccesfully() {
			if (this.equals(COMPLETED)) {
				return true;
			} else {
				return false;
			}
		}
		
		
		
		@Override
		public String toString() {
			return description;
		}
	};
	
	
	/**
	 * Encapsulates notification to make it passable to Event Dispatch Thread.
	 */
	private class TaskStateChangeNotifier implements Runnable {
		private Task.State oldState;
		private Task.State newState;

		public TaskStateChangeNotifier(Task.State oldState,	Task.State newState) {
			this.oldState = oldState;
			this.newState = newState;
		}

		public void run() {
			for (TaskEventListener listener: listeners) {
				listener.onStateChange(Task.this, oldState, newState);	
			}
		}
	}

	
	
	private OperationRecord operationRecord;
	private String id; 
	private State state = State.NEW;
	private String stateDetail = "";
	private int completionPercentage = -1;
	private long startTime;
	private long endTime;
	private String errorMessage;
	private String screenOutput;
	private String sourceCode;
	private ArrayList<DataBean> outputs = new ArrayList<>();
	private boolean hasBeenRetried = false;
	private boolean hidden = false;
	
	private List<TaskEventListener> listeners = new LinkedList<TaskEventListener>();
	private boolean isLocal;	
	
	public Task(OperationRecord operation) {
		this.operationRecord = operation;
		this.id = generateId();
//		this.isLocal = operation.getDefinition().isLocal();
	}

	public String getOperationID() {
		return operationRecord.getNameID().getID();
	}
	
	public String getName() {
		return operationRecord.getNameID().getDisplayName();
	}
	
	public String getFullName() {
		return operationRecord.getCategoryName() + " / " + operationRecord.getNameID().getDisplayName();
	}
	
	public Iterable<DataBean> getInputs() {
		return operationRecord.getInputDataBeans();
	}
	
	
	public List<String> getParameters() throws TaskException, MicroarrayException {
		List<String> parameterStrings;
		parameterStrings = new LinkedList<String>();
		for (ParameterRecord parameter: operationRecord.getParameters()) {
			parameterStrings.add(parameter.getValue());
		}
		return parameterStrings;
	}

	public int getInputCount() {
		return operationRecord.getInputs().size();
	}
	
	public void addOutput(DataBean bean) {
		this.outputs.add(bean);
	}
	
	
	/**
	 * Set the state of the Task. Also clear stateDetail field.
	 * 
	 * Listeners are notified of the state change.
	 * 
	 * Note: Only Task and TaskExecutor classes should use this method.
	 * 
	 * 
	 * 
	 * @param state
	 */
	public synchronized void  setState(State newState) {
		State oldState = this.state;
		this.state = newState;
		this.stateDetail = "";
		
		// register change event to be invoked later
		TaskStateChangeNotifier changeNotifier = new TaskStateChangeNotifier(oldState, newState);
		SwingUtilities.invokeLater(changeNotifier);
	}


	public synchronized State getState() {
		return state; 
	}

	
	
	public void setErrorMessage(String message) {
		this.errorMessage = message;
	}

	public String getErrorMessage() {
		return errorMessage;
	}

	public Iterable<DataBean> getOutputs() {
		return outputs;
	}

	public String getScreenOutput() {
		return screenOutput;
	}

	public void setScreenOutput(String screenOutput) {
		this.screenOutput = screenOutput;
	}

	/**
	 * @return state details or empty string if details not set
	 */
	public String getStateDetail() {
		if (stateDetail != null) { 
			return stateDetail;
		} else {
			return "";
		}
	}

	/**
	 * @see #getStateDetail()
	 */
	public void setStateDetail(String stateDetail) {
		this.stateDetail = stateDetail;
	}

	public String getId() {
		return id;
	}

	public void setId(String id) {
		this.id = id;
	}
	
	public void changeId() {
		this.id = generateId();
	}

	public boolean hasBeenRetried() {
		return hasBeenRetried;
	}

	public void setHasBeenRetried(boolean hasBeenRetried) {
		this.hasBeenRetried = hasBeenRetried;
	}

	public long getStartTime() {
		return startTime;
	}

	public void setStartTime(long startTime) {
		this.startTime = startTime;
	}

	public long getEndTime() {
		return endTime;
	}

	public void setEndTime(long endTime) {
		this.endTime = endTime;
	}

	/**
	 * Execution time is the time between passing a message corresponding 
	 * to this task to JMS layer and receiving+parsing a reply message.
	 */
	public long getExecutionTime() {
		if (endTime > startTime) {
			return endTime - startTime;
		} else if (startTime > 0) {
			return System.currentTimeMillis() - startTime;
		} else {
			return 0;
		}
	}

	private String generateId() {
		return UUID.randomUUID().toString();
	}

	public boolean isHidden() {
		return hidden;
	}

	public void addTaskEventListener(TaskEventListener listener) {
		listeners.add(listener);
	}


	public void setCompletionPercentage(int completionPercentage) {
		this.completionPercentage = completionPercentage;
		
	}

	public int getCompletionPercentage() {
		return completionPercentage;
	}

	public void setSourceCode(String sourceCode) {
		this.sourceCode = sourceCode;
	}

	public String getSourceCode() {
		return this.sourceCode;
	}

	public boolean isLocal() {
		return isLocal;
	}

	public OperationRecord getOperationRecord() {
		return operationRecord;
	}
}
