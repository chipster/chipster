package fi.csc.microarray.analyser;

public class JobCancelledException extends AnalysisException {

	public JobCancelledException(Exception cause) {
		super(cause);
	}

	public JobCancelledException(String message) {
		super(message);
	}

	public JobCancelledException() {
		super();
	}
}
