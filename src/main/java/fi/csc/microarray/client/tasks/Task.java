/*
 * Created on Feb 10, 2005
 *
 */
package fi.csc.microarray.client.tasks;

import java.util.ArrayList;
import java.util.Collection;
import java.util.Date;
import java.util.LinkedList;
import java.util.List;
import java.util.UUID;

import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.client.operation.OperationRecord.InputRecord;
import fi.csc.microarray.client.operation.OperationRecord.ParameterRecord;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.ThreadUtils;

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
	private Date startTime;
	private Date endTime;
	private String errorMessage;
	private String screenOutput;
	private ArrayList<DataBean> outputs = new ArrayList<>();
	private boolean hidden = false;
	
	private List<TaskEventListener> listeners = new LinkedList<TaskEventListener>();
	private boolean isLocal;	
	
	public Task(OperationRecord operationRecord, boolean local) {
		this.operationRecord = operationRecord;
		this.id = generateId();
		this.isLocal = local;
	}

	public Task(OperationRecord operationRecord, String jobId, Date startTime, Date endTime, boolean local) {
		this.operationRecord = operationRecord;
		this.id = jobId;
		this.isLocal = local;
		this.startTime = startTime;
		this.endTime = endTime;
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
	
	public Collection<InputRecord> getInputRecords() {
		return operationRecord.getInputRecords();
	}
	
	public Iterable<DataBean> getInputDataBeans() {
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
		return operationRecord.getInputRecords().size();
	}
	
	public void addOutput(DataBean bean) {
		this.outputs.add(bean);
	}
	
	
	/**
	 * Set the state of the Task. Also clear stateDetail field.
	 * 
	 * Note: Only Task and TaskExecutor classes should use this method.
	 * 
	 * Note2: Call notifyTaskStateChangeListener() after this, but outside any
	 * synchronized(Task) blocks. 
	 * 
	 * @param state
	 */
	public synchronized void  setState(State newState) {
		this.state = newState;
		this.stateDetail = "";
	}
	
	public void notifyTaskStateChangeListener(State oldState, State newState) {
		/*
		 * Notify listener
		 * 
		 * Use invokeAndWait instead of invokeLater. In CLI, we must know when
		 * the task is completed and the next operation (usually session
		 * saving) can be started.
		 */		
		TaskStateChangeNotifier changeNotifier = new TaskStateChangeNotifier(oldState, newState);
		
		ThreadUtils.runInEDT(changeNotifier);		
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

	public List<DataBean> getOutputs() {
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

	public Date getStartTime() {
		return startTime;
	}

	public void setStartTime(Date startTime) {
		this.startTime = startTime;
		if (this.operationRecord != null) {
			this.operationRecord.setStartTime(startTime);
		}
	}

	public Date getEndTime() {
		return endTime;
	}

	public void setEndTime(Date endTime) {
		this.endTime = endTime;
		if (this.operationRecord != null) {
			this.operationRecord.setEndTime(endTime);
		}
	}

	/**
	 * Execution time is the time between passing a message corresponding 
	 * to this task to JMS layer and receiving+parsing a reply message.
	 * @return current time - start time if not finished, -1 if finished but start time or end time not available
	 */
	public long getExecutionTime() {
		// finished
		if (this.state.isFinished()) {
			if (startTime != null && endTime != null) {
				return endTime.getTime() - startTime.getTime();
			}
		} 
		
		// still running
		else {
			if (startTime != null) {
				return System.currentTimeMillis() - startTime.getTime();	
			}
		}

		// something missing
		return -1;
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

	public boolean isLocal() {
		return isLocal;
	}

	public OperationRecord getOperationRecord() {
		return operationRecord;
	}
}
