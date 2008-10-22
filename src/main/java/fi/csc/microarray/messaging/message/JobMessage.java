/*
 * Created on Feb 11, 2005
 *
 *
 */
package fi.csc.microarray.messaging.message;

import java.util.List;
import javax.jms.JMSException;
import javax.jms.MapMessage;

import org.apache.log4j.Logger;


/**
 * For sending jobs to back-end components.
 * 
 * @author hupponen, akallio
 *
 */
public class JobMessage extends PayloadMessage {

	private static final Logger logger = Logger
	.getLogger(JobMessage.class);

	
	private static final String KEY_JOB_ID = "jobID";
	private static final String KEY_ANALYSIS_ID = "analysisID";
	
	private String analysisId;
	private String jobId;
	

	/**
	 * For reflection compatibility (newInstance). DO NOT REMOVE!
	 */
	public JobMessage() {
		super();
	}
	
	public JobMessage(String jobId, String analysisId, List<String> parameters) {
		super(parameters);
		this.jobId = jobId;
		this.analysisId = analysisId;
	}

	@Override
	public void unmarshal(MapMessage from) throws JMSException {
		super.unmarshal(from);

		// load ids
		this.jobId = from.getString(KEY_JOB_ID);
		this.analysisId = from.getString(KEY_ANALYSIS_ID);
		logger.debug("Unmarshalled " + KEY_JOB_ID + " : " + jobId);
		logger.debug("Unmarshalled " + KEY_ANALYSIS_ID + " : " + analysisId);
	
	}
	
	/**
	 * Construct a MapMessage that can be used to create a new JobMessage.
	 */
	@Override
	public void marshal(MapMessage mapMessage) throws JMSException {
		super.marshal(mapMessage);
		
		logger.debug("Marshalling: " + KEY_JOB_ID + " : " + this.jobId);
		logger.debug("Marshalling: " + KEY_ANALYSIS_ID + " : " + this.analysisId);
		
		// add ids
		mapMessage.setString(KEY_JOB_ID, this.jobId);
		mapMessage.setString(KEY_ANALYSIS_ID, this.analysisId);
	}
	
	/**
	 * Returns identifier of the requested job.
	 */
	public String getAnalysisId() {
		return this.analysisId;
	}

	/**
	 * @see #getAnalysisId()
	 */
	public void setAnalysisId(String id) {
		this.analysisId = id;
	}

	public String getJobId() {
		return jobId;
	}

	public void setJobId(String jobId) {
		this.jobId = jobId;
	}


}
	

