/*
 * Created on Feb 24, 2005
 *
 */
package fi.csc.microarray.comp;

import fi.csc.microarray.exception.MicroarrayException;

/**
 * Problem when running analysis jobs.
 * 
 * @author hupponen
 */
public class CompException extends MicroarrayException {

	public CompException(String message) {
		super(message);
	}

	public CompException(Exception cause) {
		super(cause);
	}

	public CompException() {
		super();
	}
}
