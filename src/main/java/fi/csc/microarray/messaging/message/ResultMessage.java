/*
 * Created on Feb 11, 2005
 *
 *
 */
package fi.csc.microarray.messaging.message;

import javax.jms.Destination;
import javax.jms.JMSException;
import javax.jms.MapMessage;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.JobState;

/**
 * For returning job results to clients.
 * 
 * @author hupponen, akallio
 *
 */
public class ResultMessage extends PayloadMessage {
	/**
	 * Logger for this class
	 */
	@SuppressWarnings("unused")
	private static final Logger logger = Logger.getLogger(ResultMessage.class);
		
	private static final String KEY_STATE = "exitState";
	private static final String KEY_JOB_ID = "jobId";
	private static final String KEY_STATE_DETAIL = "stateDetail";
	private static final String KEY_ERROR_MESSAGE = "errorMessage";
	private static final String KEY_OUTPUT_TEXT = "outputText";
	
	private String jobId;
	private JobState state;
	private String stateDetail;
	private String errorMessage;
	private String outputText;
	

	public ResultMessage(String jobId, JobState state, String stateDetail, String errorMessage,
			String outputText, 
			Destination replyTo) {
        super(null); // null because we don't have parameters
        super.setReplyTo(replyTo);
		this.jobId = jobId;
        this.state = state;
		this.stateDetail = stateDetail;
		this.errorMessage = errorMessage;
		this.outputText = outputText;
	}
	
	public ResultMessage() {
		super();
	}
	
	public void unmarshal(MapMessage from) throws JMSException {
		super.unmarshal(from);
	
		this.jobId = from.getString(KEY_JOB_ID);
		this.state = JobState.valueOf(from.getString(KEY_STATE));
		this.stateDetail = from.getString(KEY_STATE_DETAIL);
		this.errorMessage = from.getString(KEY_ERROR_MESSAGE);
		this.outputText = from.getString(KEY_OUTPUT_TEXT);
	}

	public void marshal(MapMessage mapMessage) throws JMSException {
		super.marshal(mapMessage);
		
		mapMessage.setString(KEY_JOB_ID, this.jobId);		
		mapMessage.setString(KEY_STATE, this.state.name());
		mapMessage.setString(KEY_STATE_DETAIL, this.stateDetail);
		mapMessage.setString(KEY_ERROR_MESSAGE, this.errorMessage);
		mapMessage.setString(KEY_OUTPUT_TEXT, this.outputText);
	}
	
	/**
	 * Return error message in case of failed job execution. 
	 */
	public String getErrorMessage() {
		return errorMessage;
	}

	/**
	 * @see #getErrorMessage()
	 */
	public void setErrorMessage(String errorMessage) {
		this.errorMessage = errorMessage;
	}

	/**
	 * Return the exit state of the job.
	 */
	public JobState getState() {
		return state;
	}

	/**
	 * @see #getState()
	 */
	public void setState(JobState exitState) {
		this.state = exitState;
	}

	/**
	 * Returns the text output (sysout) of the job. 
	 */
	public String getOutputText() {
		return outputText;
	}

	/**
	 * @see #getOutputText()
	 */
	public void setOutputText(String output) {
		this.outputText = output;
	}

	public String getStateDetail() {
		return stateDetail;
	}

	public void setStateDetail(String stateDetail) {
		this.stateDetail = stateDetail;
	}

	public String getJobId() {
		return jobId;
	}

	public void setJobId(String jobId) {
		this.jobId = jobId;
	}
}
	

