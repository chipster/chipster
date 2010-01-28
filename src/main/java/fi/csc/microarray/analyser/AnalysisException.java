/*
 * Created on Feb 24, 2005
 *
 */
package fi.csc.microarray.analyser;

import fi.csc.microarray.exception.MicroarrayException;

/**
 * Problem when running analysis jobs.
 * 
 * @author hupponen
 */
public class AnalysisException extends MicroarrayException {

	public AnalysisException(String message) {
		super(message);
	}

	public AnalysisException(Exception cause) {
		super(cause);
	}

	public AnalysisException() {
		super();
	}
}
