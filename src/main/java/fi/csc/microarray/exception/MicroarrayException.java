/*
 * Created on Mar 24, 2005
 *
 */
package fi.csc.microarray.exception;

/**
 * Wrapped exception for all exceptions caused by the Chipster system.
 * 
 * @author Aleksi Kallio
 */
public class MicroarrayException extends Exception {
	
	private Object extraInfo;
	
	public MicroarrayException() {
		super();
	}
	
	public MicroarrayException(String message) {
		super(message);
	}

	public MicroarrayException(Exception e) {
		this(e.getMessage(), e);
	}

	public MicroarrayException(String message, Exception e) {
		super(message, e);
	}

	public Object getExtraInfo() {
		return extraInfo;
	}

}
