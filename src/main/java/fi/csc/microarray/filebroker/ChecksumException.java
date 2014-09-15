package fi.csc.microarray.filebroker;


public class ChecksumException extends Exception {

	public ChecksumException(String message, ChecksumException cause) {
		super(message, cause);
	}

	public ChecksumException() {
		super("cheksum comparison failed: corrupted data");
	}
}
