package fi.csc.microarray.client.dialog;

public class DialogInfo {
	
	public enum Severity {
		INFO("Notification"),
		WARNING("Warning"),
		ERROR("Error"),
		QUESTION("Question");
		
		private String description;

		private Severity(String description) {
			this.description = description;
		}
		
		@Override
		public String toString() {
			return description;
		}
	}

	public enum Type {
		MESSAGE("Close"),
		OK_MESSAGE("OK"),
		OPTION("OK"),
		BLOCKER(null);
		
		private String buttonText;
		
		private Type(String buttonText) {
			this.buttonText = buttonText;
		}
		
		public String getButtonText() {
			return buttonText;
		}
	}
	
	private Severity severity = null;
	private String title = null;
	private String message = null;
	private String details = null;
    private Boolean showFeedback = false;
	private Type type;

	public DialogInfo(Severity severity, String title, String message, String details) {
		this(severity, title, message, details, Type.MESSAGE);
	}
	
	public DialogInfo(Severity severity, String title, String message, String details, Type type) {
		setSeverity(severity);
		setTitle(title);
		setMessage(message);
		setDetails(details);
		setType(type);
	}
		
	public Severity getSeverity() {
		return severity;
	}

	public String getMessage() {
		return message;
	}

	public String getDetails() {
		return details;
	}

	public void setSeverity(Severity severity) {
		this.severity = severity;
	}

	public void setMessage(String message) {
		this.message = message;
	}

	public void setDetails(String details) {
		this.details = details;
	}

	public String getTitle() {
		return title;
	}

	public void setTitle(String title) {
		this.title = title;
	}

	public Type getType() {
		return type;
	}

	public void setType(Type type) {
		this.type = type;
	}
	
	/**
	 * Set visibility of "Send report" button with which
	 * user can give us feedback.
	 * 
	 * @param show
	 */
	public void setFeedbackVisible(Boolean show) {
	    this.showFeedback = show;
	}
	
	/**
	 * Decide if "Send report" button is visible in this
	 * dialog.
	 * 
	 * @return true if feedback button is visible, false
	 * otherwise.
	 */
    public Boolean getFeedbackVisible() {
        return showFeedback;
    }
}
