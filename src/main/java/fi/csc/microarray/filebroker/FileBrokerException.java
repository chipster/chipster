package fi.csc.microarray.filebroker;


@SuppressWarnings("serial")
public class FileBrokerException extends Exception {

	public FileBrokerException() {
		super();
	}
	
	public FileBrokerException(String message) {
		super(message);
	}

	public FileBrokerException(Exception e) {
		super(e);
	}
}
