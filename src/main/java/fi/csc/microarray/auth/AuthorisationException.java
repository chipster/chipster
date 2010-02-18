/*
 * Created on Mar 24, 2005
 *
 */
package fi.csc.microarray.auth;

import fi.csc.microarray.exception.MicroarrayException;

/**
 * @author akallio
 *
 */
public class AuthorisationException extends MicroarrayException {
	
	public AuthorisationException(String s) {
		super(s);
	}

}
