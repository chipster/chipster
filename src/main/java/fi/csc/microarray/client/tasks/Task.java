/*
 * Created on Feb 10, 2005
 *
 */
package fi.csc.microarray.client.tasks;

import java.io.IOException;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.UUID;

import javax.swing.SwingUtilities;

import org.apache.log4j.Logger;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.exception.MicroarrayException;

// TODO simplify Task, takes only input and Operation
/**
 * @author Aleksi Kallio
 *
 */
public class Task {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(Task.class);
	

	public enum State {
		NEW("New"), 
		TRANSFERRING_INPUTS("Transferring inputs"),
		WAITING("Waiting"),
		RUNNING("Running"), 
		TRANSFERRING_OUTPUTS("Transferring outputs"),
		COMPLETED("Completed"), 
		FAILED("Failed"),
		FAILED_USER_ERROR("Failed, user error"),
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

	
	
	
	
	private Map<String, DataBean> inputs = new HashMap<String, DataBean>();
	private List<Object> parameters = new LinkedList<Object>();
	
	private State state = State.NEW;
	private String stateDetail = "";
	private int completionPercentage = -1;
	private String name;
	private long startTime;
	private long endTime;
	private String errorMessage;
	private String screenOutput;
	private Map<String, DataBean> outputs = new HashMap<String, DataBean>();
	private String id; 
	private boolean hasBeenRetried = false;
	private boolean hidden = false;
	
	private List<TaskEventListener> listeners = new LinkedList<TaskEventListener>();	
	
	public Task(String name) {
		this.name = name;
		this.id = generateId();
	}

	public Task(String name, boolean hidden) {
		this(name);
		this.hidden = hidden;
	}

	/**
	 * @return Returns the name.
	 */
	public String getName() {
		return name;
	}
	
	public String getNamePrettyPrinted() {
		return name.replaceAll("\"", "").replaceAll("/", " / ");
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
	
	/**
	 * A generic input method. Infers correct addInput-method to call by using 
	 * the type of input. Please use type-specific addInput-methods where 
	 * possible. 
	 * 
	 * @throws IllegalArgumentException if parameter input is of unsupported type
	 * @param name
	 * @param input
	 */
	public void addParameter(String name, Object input) {
		if (input instanceof Float) {
			addParameter(name, (Float)input);
		} else if (input instanceof Integer) {
			addParameter(name, (Integer)input);
		} else if (input instanceof String) {
			addParameter(name, (String)input);
		} else {
			throw new IllegalArgumentException("unsupported input type: " + input.getClass().getSimpleName());
		}
	}
	
	public void addParameter(String name, Integer input) {
		addParameter(name, input.toString()); // we handle Integer internally as a String
	}
	
	public void addParameter(String name, Float input) {
		addParameter(name, input.toString()); // we handle Float internally as a String
	}
	
	public void addParameter(String name, String input) {
		// we don't actually need the name now, order is enough
		logger.debug("added parameter " + name + " -> " + input);
		parameters.add(input);
	}
	
	public DataBean getInput(String name) {
		return inputs.get(name);
	}
	
	
	public List<String> getParameters() throws TaskException, MicroarrayException {
		List<String> inputArray = new LinkedList<String>();
		for (Object input : parameters) {
			if (input instanceof String) {
				inputArray.add((String)input);
			} else {
				throw new TaskException("do not know how to encode" + input.getClass().getName());
			}
		}
		return inputArray;
	}
	
	public void addOutput(String name, DataBean output) throws IOException, MicroarrayException {
		outputs.put(name, output);
	}
	
	public DataBean getOutput(String name) {
		return outputs.get(name);
	}
	
	public void addInput(String name, DataBean input) {
		this.inputs.put(name, input);
	}

	public Iterable<String> inputNames() {
		return inputs.keySet();		
	}
	
	public void setErrorMessage(String message) {
		this.errorMessage = message;
	}

	public String getErrorMessage() {
		return errorMessage;
	}

	public Iterable<String> outputNames() {
		return outputs.keySet();
	}

	public Iterable<DataBean> outputs() {
		return outputs.values();
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

	public int getInputCount() {
		return inputs.size();
	}

	public void setCompletionPercentage(int completionPercentage) {
		this.completionPercentage = completionPercentage;
		
	}

	public int getCompletionPercentage() {
		return completionPercentage;
	}
}
