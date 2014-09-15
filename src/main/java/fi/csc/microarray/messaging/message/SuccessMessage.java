package fi.csc.microarray.messaging.message;

import javax.jms.JMSException;
import javax.jms.MapMessage;

import org.apache.log4j.Logger;

import fi.csc.microarray.util.Exceptions;

/** 
 * @author hupponen
 *
 */
public class SuccessMessage extends ChipsterMessage {
	/**
	 * Logger for this class
	 */
	@SuppressWarnings("unused")
	private static final Logger logger = Logger.getLogger(SuccessMessage.class);
		
	private static final String KEY_SUCCESS = "success";
	private static final String KEY_ERROR_MESSAGE = "errorMessage";
	private static final String KEY_DETAILS = "details";
	private static final String KEY_EXCEPTION = "exception";
	
	private boolean success;
	private String errorMessage;
	private String details;
	private Exception exception;
	private String exceptionString;

	public SuccessMessage(boolean success) {
		super();
		this.success = success;
	}
	
	public SuccessMessage(boolean success, String errorMessage, String details, Exception exception) {
		super();
		this.success = success;
		this.errorMessage = errorMessage;
		this.details = details;
		this.exception = exception;
	}

	public SuccessMessage(boolean success, String errorMessage, Exception exception) {
		super();
		this.success = success;
		this.errorMessage = errorMessage;
		this.exception = exception;
	}

	public SuccessMessage(boolean success, String errorMessage, String details) {
		super();
		this.success = success;
		this.errorMessage = errorMessage;
		this.details = details;
	}
	
	public SuccessMessage(boolean success, String errorMessage) {
		super();
		this.success = success;
		this.errorMessage = errorMessage;
	}

	public SuccessMessage(boolean success, Exception exception) {
		super();
		this.success = success;
		this.exception = exception;
	}
	
	public SuccessMessage() {
		super();
	}
	
	public void unmarshal(MapMessage from) throws JMSException {
		super.unmarshal(from);
		this.success = from.getBoolean(KEY_SUCCESS);
		this.errorMessage = from.getString(KEY_ERROR_MESSAGE);
		this.details = from.getString(KEY_DETAILS);
		this.exceptionString = from.getString(KEY_EXCEPTION);
	}

	public void marshal(MapMessage mapMessage) throws JMSException {
		super.marshal(mapMessage);
		mapMessage.setBoolean(KEY_SUCCESS, success);
		mapMessage.setString(KEY_ERROR_MESSAGE, this.errorMessage);
		mapMessage.setString(KEY_DETAILS, this.details);
		if (this.exception != null) {
			mapMessage.setString(KEY_EXCEPTION, this.exception.getMessage() +  Exceptions.getStackTrace(this.exception));
		}
	}
	
	
	public boolean success() {
		return success;
	}
	
	public void setSuccess(boolean success) {
		this.success = success;
	}

	public String getErrorMessage() {
		return errorMessage;
	}

	public void setErrorMessage(String errorMessage) {
		this.errorMessage = errorMessage;
	}

	public String getDetails() {
		return details;
	}

	public void setDetails(String details) {
		this.details = details;
	}

	public void setException(Exception e) {
		this.exception = e;
	}
	
	public String getExceptionString() {
		return this.exceptionString;
	}

	public String toString() {
		return 
		
		"errorMessage: " + errorMessage + "\n" + 
		"details: " + details + "\n" +
		"exception: " + exceptionString + "\n";
	}

}
	

