/*
 * Created on Feb 11, 2005
 *
 *
 */
package fi.csc.microarray.messaging.message;

import java.util.Iterator;
import java.util.List;

import javax.jms.JMSException;
import javax.jms.MapMessage;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.ToolDescription;
import fi.csc.microarray.analyser.ToolDescription.ParameterDescription;


/**
 * For sending jobs to back-end components.
 * 
 * @author Taavi Hupponen, Aleksi Kallio
 *
 */
public class JobMessage extends PayloadMessage {

	public static interface ParameterSecurityPolicy {
		/**
		 * Checks that given value is valid from a security point of view. Analysis jobs
		 * implement this to provide context dependent checking. Typically validity depends
		 * on the type of value (numeric, text...), so ParameterDescription is also passed.
		 * 
		 * @return true iff is valid
		 */
		public boolean isValueValid(String value, ParameterDescription parameterDescription);
	}
	
	@SuppressWarnings("serial")
	public static class ParameterValidityException extends Exception {
		
		public ParameterValidityException(String msg) {
			super(msg);
		}
	}
	
	private static final Logger logger = Logger.getLogger(JobMessage.class);
	
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

	/**
	 * Gets parameters in the order they were inserted.
	 * Parameters are given by the user and hence  
	 * safety policy is required to get access to them.
	 * 
	 * @param securityPolicy security policy to check parameters against, cannot be null
	 * @param description description of the analysis operation, cannot be null
	 * 
	 * @throws ParameterValidityException if some parameter value fails check by security policy 
	 */
	public List<String> getParameters(ParameterSecurityPolicy securityPolicy, ToolDescription description) throws ParameterValidityException {
		
		// Do argument checking first
		if (securityPolicy == null) {
			throw new IllegalArgumentException("security policy cannot be null");
		}
		if (description == null) {
			throw new IllegalArgumentException("analysis description cannot be null");
		}

		// Count parameter descriptions
		int parameterDescriptionCount = 0;
		for (Iterator<ParameterDescription> iterator = description.getParameters().iterator(); iterator.hasNext(); iterator.next()) {
			parameterDescriptionCount++;
		}

		// Get the actual values
		List<String> parameters = super.getParameters();

		// Check that description and values match
		if (parameterDescriptionCount != parameters.size()) {
			throw new IllegalArgumentException("number of parameter descriptions does not match the number of parameter values");
		}

		// Validate parameters
		Iterator<ParameterDescription> descriptionIterator = description.getParameters().iterator();
		for (String parameter : parameters) {
			ParameterDescription parameterDescription = descriptionIterator.next();
			if (!securityPolicy.isValueValid(parameter, parameterDescription)) {
				throw new ParameterValidityException("illegal value for parameter " + parameterDescription.getName() + ": " + parameter);
			}
		}
		
		// Everything was ok, return the parameters
		return parameters;
	}
}
	

