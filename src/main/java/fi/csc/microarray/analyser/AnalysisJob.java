package fi.csc.microarray.analyser;

import java.util.Date;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.messaging.message.ResultMessage;

/**
 * Interface to analysis jobs. Implementations do the actual analysis
 * operation.
 * 
 * @author akallio
 */
public abstract class AnalysisJob implements Runnable {

	private static final Logger logger = Logger.getLogger(AnalysisJob.class);

	
	protected JobMessage inputMessage;
	protected ResultCallback resultHandler;
	protected AnalysisDescription analysis;

	private Date receiveTime;
	private Date scheduleTime;
	private Date executionStartTime;
	private Date executionEndTime;
	
	
	
	private boolean constructed = false;

	private JobState state;
	private String stateDetail;
	private boolean toBeCanceled = false;
	protected ResultMessage outputMessage;
	
	
	
	
	
	public AnalysisJob() {
		this.state = JobState.NEW;
		this.stateDetail = "New job created.";
	}
	

	public void construct(JobMessage inputMessage, AnalysisDescription analysis, ResultCallback resultHandler) {
		this.constructed = true;
		this.analysis = analysis;
		this.inputMessage = inputMessage;
		this.resultHandler = resultHandler;
		
		// initialize result message
		outputMessage = new ResultMessage();
		outputMessage.setJobId(this.getId());
		outputMessage.setReplyTo(inputMessage.getReplyTo());
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
		
		
			preExecute();
			this.setExecutionStartTime(new Date());
			execute();
			this.setExecutionEndTime(new Date());
			postExecute();
			
			// if we get here, the job has been successful
			
			// sanity check
			if (!getState().equals(JobState.RUNNING)) {
				throw new RuntimeException("Unexpected job end state: " + getState());
			}
			
			// update state and send results
			updateState(JobState.COMPLETED, "", false);
			outputMessage.setState(this.state);
			outputMessage.setStateDetail(this.stateDetail);
			resultHandler.sendResultMessage(inputMessage, outputMessage);
			
		} 
		// job cancelled, do nothing
		catch (JobCancelledException jce) {
			this.setExecutionEndTime(new Date());
			logger.debug("Job cancelled: " + this.getId());
		} 
		
		// send retry request for the job
		// some of the inputs were not available in the fileserver cache
		catch (RetryException re) {
			this.setExecutionEndTime(new Date());
			logger.info("Sending retry request for " + this.getId());
			updateState(JobState.RETRY, "", false);
			outputMessage.setState(this.state);
			outputMessage.setStateDetail("Waiting for resend.");
			resultHandler.sendResultMessage(inputMessage, outputMessage);
		}
		
		// job error
		catch (Throwable e) {
			this.setExecutionEndTime(new Date());
			String msg = analysis.getFullName() + " failed: ";
			
			// problem with R
			if (getState() == JobState.FAILED) {
				msg = msg + " R finished with error";
				logger.info(msg);
				outputMessage.setErrorMessage(msg);
			} 
			
			// problem caused by the user
			else if (getState() == JobState.FAILED_USER_ERROR) {
				msg = msg + " R finished with error";
				logger.info(msg);
				outputMessage.setErrorMessage(msg);
			} 
			
			// R timeout
			else if (getState() == JobState.TIMEOUT) {
				msg = msg + "R running over timeout";
				logger.info(msg);
				outputMessage.setErrorMessage(msg);
			}
			
			// unexpected server error
			else {
				msg = msg + "unexpected error";
				logger.error(msg, e);
				updateState(JobState.ERROR, "Unexpected server error", false);
				outputMessage.setErrorMessage(e.toString());
			}

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
	
	public synchronized void updateState(JobState newState, String newStateDetail, boolean sendNotification) {

		if (getState() == JobState.CANCELLED) {
			return;
		}

		
		// update state
		this.state = newState;
		this.stateDetail = newStateDetail;
	
		// send notification message
		if (sendNotification) {
			outputMessage.setState(this.state);
			outputMessage.setStateDetail(this.stateDetail);
			resultHandler.sendResultMessage(inputMessage, outputMessage);
		}
	}

	public synchronized void updateStateDetail(String newStateDetail, boolean sendNotification) {
		
		if (this.state.equals(JobState.CANCELLED)) {
			return;
		}
		
		// update state
		this.stateDetail = newStateDetail;
	
		// send notification message
		if (sendNotification) {
			outputMessage.setState(this.state);
			outputMessage.setStateDetail(this.stateDetail);
			resultHandler.sendResultMessage(inputMessage, outputMessage);
		}
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
			updateState(JobState.CANCELLED, "", false);
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

	
	
	protected abstract void execute() throws Exception;
	
	
	
	protected void preExecute() throws Exception {
		if (!constructed) {
			throw new IllegalStateException("you must call construct(...) first");
		}
		updateState(JobState.RUNNING, "Initialising", true);
	}
	
	protected void postExecute()  throws Exception {
		
	}

	protected void cleanUp() {
	}


	public JobMessage getInputMessage() {
		return inputMessage;
	}

	public ResultMessage getResultMessage() {
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

}