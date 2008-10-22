package fi.csc.microarray;

import fi.csc.microarray.util.Exceptions;

public class ErrorReportAsException extends MicroarrayException {

	private String title = null;
	private String message = null;
	private String details = null;
	
	
	public ErrorReportAsException(String title, String message, String details) {
		super(title);
		this.title = title;
		this.message = message;
		this.details = details;
	}


	public ErrorReportAsException(String title, String message, Exception e) {
		this(title, message, Exceptions.getStackTrace(e));
	}


	public String getTitle() {
		return title;
	}


	public String getMessage() {
		return message;
	}


	public String getDetails() {
		return details;
	}
	
	
}
