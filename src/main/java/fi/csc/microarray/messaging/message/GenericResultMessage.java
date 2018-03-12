package fi.csc.microarray.messaging.message;

import java.time.Instant;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.JobState;

/**
 * Generic result message, independent of the communication method used.
 * 
 * @author klemela
 *
 */
public class GenericResultMessage {
	/**
	 * Logger for this class
	 */
	@SuppressWarnings("unused")
	private static final Logger logger = Logger.getLogger(GenericResultMessage.class);
		
	private String jobId;
	private JobState state;
	private String stateDetail;
	private String errorMessage;
	private String outputText;
	private String sourceCode;
	private Instant startTime;
	private Instant endTime;
	
	private Map<String, String> ids = new HashMap<String, String>();
	private Map<String, String> names = new HashMap<String, String>();
	
	public GenericResultMessage(String jobId, JobState state, String stateDetail, String errorMessage,
			String outputText) {
		
		this.jobId = jobId;
        this.state = state;
		this.stateDetail = stateDetail;
		this.errorMessage = errorMessage;
		this.outputText = outputText;
	}
	
	public GenericResultMessage() {
		super();
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

	public void setSourceCode(String sourceCode) {
		this.sourceCode = sourceCode;
	}
	
	public String getSourceCode() {
		return this.sourceCode;
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
	
	public void addDataset(String outputName, String id, String name) {
		ids.put(outputName, id);
		names.put(outputName, name);
	}
	
	public Set<String> getOutputNames() {
		HashSet<String> keys = new HashSet<>();
		keys.addAll(ids.keySet());
		keys.addAll(names.keySet());
		return keys;
	}
	
	public String getDatasetId(String outputName) {
		return ids.get(outputName);
	}

	public String getDatasetName(String outputName) {
		return names.get(outputName);
	}

	public Instant getStartTime() {
		return startTime;
	}

	public void setStartTime(Instant startTime) {
		this.startTime = startTime;
	}

	public Instant getEndTime() {
		return endTime;
	}

	public void setEndTime(Instant endTime) {
		this.endTime = endTime;
	}


}
	

