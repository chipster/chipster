/*
 * Created on Feb 11, 2005
 *
 *
 */
package fi.csc.microarray.messaging.message;

import java.text.DateFormat;
import java.text.ParseException;
import java.util.Date;
import java.util.HashMap;

import javax.jms.JMSException;
import javax.jms.MapMessage;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.JobState;

/** 
 * @author hupponen
 *
 */
public class JobLogMessage extends ChipsterMessage {
	/**
	 * Logger for this class
	 */
	@SuppressWarnings("unused")
	private static final Logger logger = Logger.getLogger(JobLogMessage.class);
		
	private static final String KEY_OPERATION = "operation";
	private static final String KEY_STATE = "exitState";
	private static final String KEY_STATE_DETAIL = "stateDetail";
	private static final String KEY_JOB_ID = "jobId";
	private static final String KEY_START_TIME = "startTime";
	private static final String KEY_END_TIME = "endTime";
	private static final String KEY_ERROR_MESSAGE = "errorMessage";
	private static final String KEY_OUTPUT_TEXT = "outputText";
	
	private static final String KEY_USERNAME = "username";
	private static final String KEY_COMP_HOST = "compHost";
	
	private String operation;
	private JobState state;
	private String stateDetail; // not logged
	private String jobId;
	private Date startTime;
	private Date endTime;
	private String errorMessage;
	private String outputText;
	private String username;
	private String compHost;

	
	public JobLogMessage(String operation, JobState state, String stateDetail, String jobId, Date startTime, Date endTime, String errorMessage, String outputText, String username, String compHost) {
		super();
		this.operation = operation;
		this.state = state;
		this.stateDetail = stateDetail;
		this.jobId = jobId;
		this.startTime = startTime;
		this.endTime = endTime;
		this.errorMessage = errorMessage;
		this.outputText = outputText;
		this.username = username;
		this.compHost = compHost;
	}
	
	public JobLogMessage() {
		super();
	}
	
	public void unmarshal(MapMessage from) throws JMSException {
		super.unmarshal(from);
	
		this.operation = from.getString(KEY_OPERATION);
		if (from.getString(KEY_STATE) != null) {
			this.state = JobState.valueOf(from.getString(KEY_STATE));
		}
		
		this.stateDetail = from.getString(KEY_STATE_DETAIL);
		
		this.jobId = from.getString(KEY_JOB_ID);
		try {
			DateFormat df = DateFormat.getDateTimeInstance();
			if (from.getString(KEY_START_TIME) != null) {
				this.startTime = df.parse(from.getString(KEY_START_TIME));
			}
			if (from.getString(KEY_END_TIME) != null) {
				this.endTime = df.parse(from.getString(KEY_END_TIME));
			}
		} catch (ParseException e) {
			throw new JMSException(e.toString());
		}
		this.errorMessage = from.getString(KEY_ERROR_MESSAGE);
		this.outputText = from.getString(KEY_OUTPUT_TEXT);
		this.username = from.getString(KEY_USERNAME);
		this.compHost = from.getString(KEY_COMP_HOST);
	}

	public void marshal(MapMessage mapMessage) throws JMSException {
		super.marshal(mapMessage);
		
		mapMessage.setString(KEY_OPERATION, this.operation);
		mapMessage.setString(KEY_STATE, this.state.name());
		mapMessage.setString(KEY_STATE_DETAIL, this.stateDetail);
		mapMessage.setString(KEY_JOB_ID, this.jobId);		
		
		DateFormat df = DateFormat.getDateTimeInstance();
		if (this.startTime != null) {
			mapMessage.setString(KEY_START_TIME, df.format(this.startTime));
		}
		if (this.endTime != null) {
			mapMessage.setString(KEY_END_TIME, df.format(this.endTime));
		}
		
		mapMessage.setString(KEY_ERROR_MESSAGE, this.errorMessage);
		mapMessage.setString(KEY_OUTPUT_TEXT, this.outputText);
		mapMessage.setString(KEY_USERNAME, this.username);
		mapMessage.setString(KEY_COMP_HOST, this.compHost);
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

	public String getJobId() {
		return jobId;
	}

	public void setJobId(String jobId) {
		this.jobId = jobId;
	}

	public String getOperation() {
		return operation;
	}

	public void setOperation(String operation) {
		this.operation = operation;
	}

	public Date getStartTime() {
		return startTime;
	}

	public void setStartTime(Date startTime) {
		this.startTime = startTime;
	}

	public Date getEndTime() {
		return endTime;
	}

	public void setEndTime(Date endTime) {
		this.endTime = endTime;
	}

	public String getUsername() {
		return username;
	}

	public void setUsername(String username) {
		this.username = username;
	}

	public String getCompHost() {
		return compHost;
	}

	public void setCompHost(String compHost) {
		this.compHost = compHost;
	}

	public String toString() {
		return 
		
		"operation: " + operation + "\n" + 
		"state: " + (state != null ? state.toString() : "") + "\n" +
		"jobId: " + jobId + "\n" +
		"start time: " + (startTime != null ? startTime.toString() : "") + "\n" +
		"end time: " + (endTime != null ? endTime.toString() : "") + "\n" +
		"username: " + username + "\n" +
		"compHost: " + compHost + "\n";
	}

	public String getStateDetail() {
		return this.stateDetail;
	}

	public HashMap<String, Object> toMap() {
		HashMap<String, Object> map = new HashMap<>();
		
		DateFormat df = DateFormat.getDateTimeInstance();
		
		map.put(KEY_OPERATION, getOperation());
		map.put(KEY_STATE, getState().name());
		map.put(KEY_STATE_DETAIL, getStateDetail());
		map.put(KEY_JOB_ID, getJobId());
		
		if (startTime != null) {
			map.put(KEY_START_TIME, df.format(getStartTime()));
		}
		
		if (endTime != null) {
			map.put(KEY_END_TIME, df.format(getEndTime()));
		}
		
		map.put(KEY_ERROR_MESSAGE, getErrorMessage());
		map.put(KEY_OUTPUT_TEXT, getOutputText());
		map.put(KEY_USERNAME, getUsername());
		map.put(KEY_COMP_HOST, getCompHost());
		
		return map;
	}

	public void fromMap(HashMap<String, String> map) throws JMSException {
		
		this.operation = map.get(KEY_OPERATION);
		if (map.get(KEY_STATE) != null) {
			this.state = JobState.valueOf(map.get(KEY_STATE));
		}
		
		this.stateDetail = map.get(KEY_STATE_DETAIL);
		
		this.jobId = map.get(KEY_JOB_ID);
		try {
			DateFormat df = DateFormat.getDateTimeInstance();
			if (map.get(KEY_START_TIME) != null) {
				this.startTime = df.parse(map.get(KEY_START_TIME));
			}
			if (map.get(KEY_END_TIME) != null) {
				this.endTime = df.parse(map.get(KEY_END_TIME));
			}
		} catch (ParseException e) {
			throw new JMSException(e.toString());
		}
		this.errorMessage = map.get(KEY_ERROR_MESSAGE);
		this.outputText = map.get(KEY_OUTPUT_TEXT);
		this.username = map.get(KEY_USERNAME);
		this.compHost = map.get(KEY_COMP_HOST);
		
	}
}
	

