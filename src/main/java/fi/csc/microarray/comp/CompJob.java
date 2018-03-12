package fi.csc.microarray.comp;

import java.time.Instant;
import java.util.Date;

import org.apache.commons.lang3.StringUtils;
import org.apache.log4j.Logger;

import fi.csc.microarray.config.DirectoryLayout;
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

	public static String SCRIPT_SUCCESSFUL_STRING = "chipster-script-finished-succesfully";
	public static final String CHIPSTER_NOTE_TOKEN = "CHIPSTER-NOTE:"; 

	private static final Logger logger = Logger.getLogger(CompJob.class);
	
	protected GenericJobMessage inputMessage;
	protected ResultCallback resultHandler;
	protected ToolDescription toolDescription;

	private Date receiveTime;
	private Date scheduleTime;
	
	private int jobTimeout;
	
	

	private boolean constructed = false;

	private JobState state;
	private String stateDetail;
	private boolean toBeCanceled = false;
	private GenericResultMessage outputMessage;
	
	public CompJob() {
		outputMessage = new GenericResultMessage();
		this.state = JobState.NEW; // updateState would check old state -> NPE
		this.jobTimeout = DirectoryLayout.getInstance().getConfiguration().getInt("comp", "job-timeout");
	}
	

	public void construct(GenericJobMessage inputMessage, ToolDescription analysis, ResultCallback resultHandler) {
		this.constructed = true;
		this.toolDescription = analysis;
		this.inputMessage = inputMessage;
		this.resultHandler = resultHandler;
		
		// initialize result message
		outputMessage.setJobId(this.getId());
	}
	
	
	/**
	 * Run the job. After execution, job should report results 
	 * through the ResultCallback interface.
	 */
	public void run() {
		try {
			this.outputMessage.setStartTime(Instant.now());
			updateState(JobState.RUNNING, "initialising");
			
			if (!constructed) {
				throw new IllegalStateException("you must call construct(...) first");
			}
			
			// before execute
			preExecute();

			// execute
			if (this.getState() == JobState.RUNNING) {
				execute();
			}

			// after execute
			if (this.getState() == JobState.RUNNING) {
				postExecute();
			}

			// successful job, failed states are set when they happen 
			if (this.getState() == JobState.RUNNING) {
				updateState(JobState.COMPLETED, "");
			}
		} 

		// job was cancelled, do nothing, state has already been set when calling cancel()
		catch (JobCancelledException jce) {
			logger.debug("job cancelled: " + this.getId());
		} 
		
		
		// something unexpected happened
		catch (Throwable e) {
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

	
	
	public synchronized void updateState(JobState newState) {
		updateState(newState, "");
	}

	
	public synchronized void updateState(JobState newState, String stateDetail) {

		// ignore if jos is cancelled already
		if (this.state == JobState.CANCELLED) {
			return;
		}
		
		// should not try to update state if already finished
		if (this.state.isFinished()) {
			logger.warn("trying to update state for already finished job, old state: " + this.state.toString() + 
					"new state: " + newState.toString());
			return;
		}
		
		// set end time if new state is finished
		if (newState.isFinished()) {
			this.outputMessage.setEndTime(Instant.now());
		}
		
		// update state
		this.state = newState;
		this.stateDetail = stateDetail;

		// send notification message
		outputMessage.setState(this.state);
		outputMessage.setStateDetail(this.stateDetail);
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
	 * Check if the job should be canceled (cancel() has been called)
	 * cancels the job by throwing the cancellation exception.
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


	public Instant getStartTime() {
		return outputMessage.getStartTime();
	}

	public Instant getEndTime() {
		return outputMessage.getEndTime();
	}

	public String getStateDetail() {
		return this.stateDetail;
	}
	
	public void setSourceCode(String source) {
		this.outputMessage.setSourceCode(source);
	}
	
	public String getErrorMessage() {
		return this.outputMessage.getErrorMessage();
	}
	
	public void setErrorMessage(String message) {
		this.outputMessage.setErrorMessage(message);
	}
	
	public void setOutputText(String output) {
		this.outputMessage.setOutputText(output);
	}
	
	public void appendOutputText(String s) {
		String currentOutput = this.outputMessage.getOutputText();
		if (currentOutput != null) {
			this.outputMessage.setOutputText(currentOutput + "\n" + s);
		} else {
			this.outputMessage.setOutputText(s);
		}
	}
	
	public ToolDescription getToolDescription() {
		return toolDescription;
	}
	
	/**
	 * @return Process object of the external process or null, if this Job doesn't have any
	 */
	public Process getProcess() {
		return null;
	}
	
	protected int getTimeout() {
		return this.jobTimeout;
	}

	protected void addOutputDataset(String outputName, String id, String name) {
		outputMessage.addDataset(outputName, id, name);
	}

	
	protected static String getErrorMessage(String screenOutput, String errorMessageToken, String removeLastLineToken) {

		// find the error token
		int errorTokenStartIndex = screenOutput.lastIndexOf(errorMessageToken);

		if (errorTokenStartIndex != -1) {
			String errorMessage = screenOutput.substring(errorTokenStartIndex);

			// remove the line that contains the error token
			errorMessage = StringUtils.substringAfter(errorMessage, errorMessageToken);

			// remove last line if contains last line to remove token
			if (removeLastLineToken != null) {
				errorMessage = StringUtils.substringBeforeLast(errorMessage, removeLastLineToken);
			}

			return errorMessage.trim();
		} else {
			return null;
		}
	}

	protected static String getChipsterNote(String errorMessage) {
		// check for chipster note
		if (errorMessage.contains(CHIPSTER_NOTE_TOKEN)) {
			return errorMessage.substring(errorMessage.indexOf(CHIPSTER_NOTE_TOKEN) + CHIPSTER_NOTE_TOKEN.length())
					.trim();
		} else {
			return null;
		}
	}

	
	
}
