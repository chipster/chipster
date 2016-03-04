package fi.csc.microarray.comp;

import java.util.Date;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.GenericJobMessage;
import fi.csc.microarray.messaging.message.GenericResultMessage;
import fi.csc.microarray.util.Exceptions;

/**
 * Interface to analysis jobs. Implementations do the actual analysis
 * operation.
 * 
 * Notes for subclassing AnalysisJob:
 * 
 * Job state
 * ---------
 * Use updateState(...) and updateStateDetailToClient(...) for changing the 
 * state of the job.
 * <p>
 * If there isn't a really good reason to do otherwise, the state of a 
 * successful job should be kept as RUNNING until the AnalysisJob sets
 * it to COMPLETED. So don't set the state as COMPLETED in a subclass.
 * <p>
 * Try to detect any error situations and set the state accordingly:
 * <p>
 * FAILED_USER_ERROR:	A situation in which the user can fix the problem
 * 						by changing inputs or parameters. You must also
 * 						set the error message to something which will
 * 						tell the user what happened and how to fix it.
 * 						The error message will be used as the header
 * 						in the client dialog.
 * <p>
 * FAILED:				The tool failed.
 * <p>
 * ERROR:				Error in the system, not in a specific tool.
 * <p>
 * Before setting the state to one of those above, set the error message
 * and output text, then return after updateState(...)
 * <p>
 * AnalysisJob will catch any Exceptions (and Throwables) and set the 
 * state as ERROR (With the exception of JobCancelledException). 
 * <p>
 * Call cancelCheck() in situations where it makes sense to acknowledge
 * that the jobs has been requested to be canceled. This will throw the
 * JobCancelledException.
 * <p>
 * Use updateStateDetailToClient(...) for reporting meaningful state 
 * detail changes to the client.
 * 
 * 
 * @author akallio, hupponen
 */
public abstract class CompJob implements Runnable {

	public static String SCRIPT_SUCCESSFUL_STRING = "script-finished-succesfully";
	public static String SCRIPT_FAILED_STRING = "script-finished-unsuccesfully";	
	public final String CHIPSTER_NOTE_TOKEN = "CHIPSTER-NOTE:"; 

	private static final Logger logger = Logger.getLogger(CompJob.class);
	
	protected GenericJobMessage inputMessage;
	protected ResultCallback resultHandler;
	protected ToolDescription toolDescription;

	private Date receiveTime;
	private Date scheduleTime;
	private Date executionStartTime;
	private Date executionEndTime;

	private boolean constructed = false;

	private JobState state;
	private String stateDetail;
	private boolean toBeCanceled = false;
	protected GenericResultMessage outputMessage;
	
	public CompJob() {
		this.state = JobState.NEW;
		this.stateDetail = "new job created.";
	}
	

	public void construct(GenericJobMessage inputMessage, ToolDescription analysis, ResultCallback resultHandler) {
		this.constructed = true;
		this.toolDescription = analysis;
		this.inputMessage = inputMessage;
		this.resultHandler = resultHandler;
		
		// initialize result message
		outputMessage = new GenericResultMessage();
		outputMessage.setJobId(this.getId());
		outputMessage.setState(this.state);
		outputMessage.setStateDetail(this.stateDetail);
	}
	
	
	/**
	 * Run the job. After execution, job should report results 
	 * through the ResultCallback interface.
	 */
	public void run() {
		try {
			if (!constructed) {
				throw new IllegalStateException("you must call construct(...) first");
			}

			// before execute
			preExecute();

			// execute
			if (this.getState() == JobState.RUNNING) {
				this.setExecutionStartTime(new Date());
				execute();
				this.setExecutionEndTime(new Date());
			}

			// after execute
			if (this.getState() == JobState.RUNNING) {
				postExecute();
			}

			// successful job
			if (this.getState() == JobState.RUNNING) {
				// update state and send results
				this.state = JobState.COMPLETED;
				this.stateDetail = "";
			}
			
			// send the result message 
			// for the unsuccessful jobs, the state has been set by the subclass
			outputMessage.setState(this.state);
			outputMessage.setStateDetail(this.stateDetail);
			resultHandler.sendResultMessage(inputMessage, outputMessage);
			
		} 
		// job cancelled, do nothing
		catch (JobCancelledException jce) {
			this.setExecutionEndTime(new Date());
			logger.debug("job cancelled: " + this.getId());
		} 
		
		
		// something unexpected happened
		catch (Throwable e) {
			this.setExecutionEndTime(new Date());
			updateState(JobState.ERROR, "running tool failed");
			outputMessage.setErrorMessage("Running tool failed.");
			outputMessage.setOutputText(Exceptions.getStackTrace(e));
			outputMessage.setState(this.state);
			outputMessage.setStateDetail(this.stateDetail);
			resultHandler.sendResultMessage(inputMessage, outputMessage);
		} 

		// clean up
		finally {
			try {
				cleanUp();
			} catch (Throwable t) {
				logger.error("Error when cleaning up.", t);
			}
			resultHandler.removeRunningJob(this);
		}
	}

	public String getId() {
		return this.inputMessage.getJobId();
	}
	
	public synchronized void updateState(JobState newState, String stateDetail) {

		// don't allow new state changes, if this is cancelled already
		if (getState() == JobState.CANCELLED) {
			return;
		}
		
		// update state
		this.state = newState;
		this.stateDetail = stateDetail;
	}

	/**
	 * Updates the state detail and also sends the state and the detail 
	 * to the client.
	 * 
	 * @param newStateDetail
	 */
	public synchronized void updateStateDetailToClient(String newStateDetail) {

		// job may continue for some time before it checks if it's cancelled
		if (getState() == JobState.CANCELLED) {
			return;
		}

		// update state
		this.stateDetail = newStateDetail;

		// send notification message
		outputMessage.setState(this.state);
		outputMessage.setStateDetail(this.stateDetail);
		resultHandler.sendResultMessage(inputMessage, outputMessage);
	}


	public synchronized void updateStateToClient() {
		this.updateStateToClient(this.state, this.stateDetail, true);
	}
	
	public synchronized void updateStateToClient(JobState newState, String stateDetail) {
		this.updateStateToClient(newState, stateDetail, false);
	}

	public synchronized void updateStateToClient(JobState newState, String stateDetail, boolean isHeartbeat) {
		
		// job may continue for some time before it checks if it's cancelled
		if (getState() == JobState.CANCELLED) {
			return;
		}
		
		// update state
		this.state = newState;
		this.stateDetail = stateDetail;

		// send notification message
		outputMessage.setState(this.state);
		outputMessage.setStateDetail(this.stateDetail);
		outputMessage.setHeartbeat(isHeartbeat);
		resultHandler.sendResultMessage(inputMessage, outputMessage);
	}	
	
	public JobState getState() {
		return this.state;
	}
	
	
	/**
	 * Request the job to be canceled. The job is not canceled immediately. Instead it is 
	 * flagged to be canceled and will cancel as soon as possible.
	 * 
	 */
	public void cancel() {
		logger.debug("Canceling job " + getId());
		this.toBeCanceled = true;
		updateState(JobState.CANCELLED, "");
		cancelRequested();
	}
	
	
	/**
	 * Check if the job should be canceled (cancel() has been called), and if needed
	 * updates the state and cancels the job by throwing the cancellation exception.
	 * 
	 * Should be called by subclasses in situations where it is safe to cancel the job.
	 * 
	 * 
	 * @throws JobCancelledException
	 */
	protected void cancelCheck() throws JobCancelledException {
		if (toBeCanceled) {			
			throw new JobCancelledException();
		}
	}

	/**
	 * Will be called when the job is canceled using the cancel() method.
	 * 
	 * Subclasses should override this method with actions that should be taken
	 * immediately when cancel is requested.
	 *
	 */
	protected abstract void cancelRequested();

	
	
	protected abstract void execute() throws JobCancelledException;
	
	
	
	protected void preExecute() throws JobCancelledException {
		if (!constructed) {
			throw new IllegalStateException("you must call construct(...) first");
		}
		updateStateToClient(JobState.RUNNING, "initialising");
	}
	
	protected void postExecute()  throws JobCancelledException {
		
	}

	protected void cleanUp() {
	}


	public GenericJobMessage getInputMessage() {
		return inputMessage;
	}

	public GenericResultMessage getResultMessage() {
		return outputMessage;
	}

	public Date getReceiveTime() {
		return receiveTime;
	}


	public void setReceiveTime(Date receiveTime) {
		this.receiveTime = receiveTime;
	}


	public Date getScheduleTime() {
		return scheduleTime;
	}


	public void setScheduleTime(Date scheduleTime) {
		this.scheduleTime = scheduleTime;
	}


	public Date getExecutionStartTime() {
		return executionStartTime;
	}


	public void setExecutionStartTime(Date executionStartTime) {
		this.executionStartTime = executionStartTime;
	}


	public Date getExecutionEndTime() {
		return executionEndTime;
	}


	public void setExecutionEndTime(Date executionEndTime) {
		this.executionEndTime = executionEndTime;
	}


	public String getStateDetail() {
		return this.stateDetail;
	}
	
	public ToolDescription getToolDescription() {
		return toolDescription;
	}
}
