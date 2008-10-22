package fi.csc.microarray.auth;

import java.io.IOException;
import java.io.Reader;
import java.util.Map;

import javax.security.auth.Subject;
import javax.security.auth.callback.Callback;
import javax.security.auth.callback.CallbackHandler;
import javax.security.auth.callback.NameCallback;
import javax.security.auth.callback.PasswordCallback;
import javax.security.auth.callback.UnsupportedCallbackException;
import javax.security.auth.login.FailedLoginException;
import javax.security.auth.login.LoginException;
import javax.security.auth.spi.LoginModule;

/**
 * 
 * @author Taavi Hupponen
 *
 */
public abstract class LoginModuleBase implements LoginModule {

	// initial state
	protected Subject subject;
	protected CallbackHandler callbackHandler;
	protected Map sharedState;
	protected Map options;
	
	// the authentication status
	private boolean succeeded = false;
	private boolean commitSucceeded = false;

	// stored info  
	private String username;

	
	
	
	public void initialize(Subject subject, CallbackHandler callbackHandler,
			Map<String, ?> sharedState, Map<String, ?> options) {
		this.subject = subject;
		this.callbackHandler = callbackHandler;
		this.sharedState = sharedState;
		this.options = options;
	}


	/**
	 *
	 * @return true in all cases since this <code>LoginModule</code>
	 *		should not be ignored.
	 *
	 * @exception FailedLoginException if the authentication fails. <p>
	 *
	 * @exception LoginException if this <code>LoginModule</code>
	 *		is unable to perform the authentication, for example if
	 *		username or password callback return null.
	 */
	public boolean login() throws LoginException {
	
		// get the username and password
		if (callbackHandler == null) { 
			throw new LoginException("No CallbackHandler.");
		}
	
		// initialise and invoke callbacks to get the username and password
		Callback[] callbacks = new Callback[2];
		NameCallback nameCallback = new NameCallback("username: ");
		PasswordCallback passwordCallback = new PasswordCallback("password: ", false);
		callbacks[0] = nameCallback;
		callbacks[1] = passwordCallback;
	    try {
			callbackHandler.handle(callbacks);
		} catch (IOException ioe) {
			throw new LoginException("Could not get username and password due to IO error.");
		} catch (UnsupportedCallbackException e) {
			throw new LoginException("Could not get username and password.");
		}
	    
	    // get the username
	    username = nameCallback.getName();
	    if (username == null) {
	    	throw new LoginException("Username is null.");
	    }
	    
	    // get the password 
	    char[] tempPassword = passwordCallback.getPassword();
	    if (tempPassword == null) {
	    	throw new LoginException("Password is null.");
	    }
	    char[] password = new char[tempPassword.length];
	    System.arraycopy(tempPassword, 0,
			password, 0, tempPassword.length);
	    passwordCallback.clearPassword();
		
		// verify username and password
	    boolean authSuccessful;
	    try {
			authSuccessful = authenticate(username, password);
		} catch (IOException e1) {
			throw new LoginException("Could not verify username and password for " + username);
		} finally {
			// clear our copy of password
			for (int i = 0; i < password.length; i++) {
				password[i] = ' ';
			}
			password = null;
		}
	
		// authentication successful
	    if (authSuccessful) {
	    	this.succeeded = true;
	    	return true;
	    } 
	    
	    // authentication failed
	    else {
	    	// clean the state
	    	username = null;
	    
	    	// finish
	    	throw new FailedLoginException("Login failed.");
	    }	    
	}


	/**
	 * This method is called if the LoginContext's overall authentication succeeded
	 * (the relevant REQUIRED, REQUISITE, SUFFICIENT and OPTIONAL LoginModules
	 * succeeded).
	 *
	 * If this LoginModule's own authentication attempt succeeded 
	 * (checked by retrieving the private state saved by the
	 * <code>login</code> method), then this method associates the Subject with the
	 * needed Principals
	 * 
	 * If this LoginModule's own authentication attempted failed, then this method removes
	 * any state that was originally saved.
	 *
	 * @exception LoginException if the commit fails.
	 *
	 * @return true if this LoginModule's own login and commit
	 *		attempts succeeded, or false otherwise.
	 */
	public boolean commit() throws LoginException {
		if (succeeded) {
			// add principals to the subject
		} 
		
		// clean out the state
		username = null;
	
		// store commit state
		commitSucceeded = true;
		
		if (succeeded && !commitSucceeded) {
			throw new LoginException("Commit failed.");
		}
		
		return succeeded && commitSucceeded;
	}


	public boolean abort() throws LoginException {
		// clean the state
		username = null;
		
		// our login successful
		if (succeeded) {
			return true;
		} else {
			return false;
		}

	}

	
	public boolean logout() throws LoginException {
		// we assign no pricipals atm, so nothing to remove from subject
		return true;
	}
	
	protected abstract boolean authenticate(String username, char[] password) throws IOException;
	
	protected static void skipToLineEnd(Reader reader) throws IOException {
		int next;
		for (next = reader.read(); next != -1; next = reader.read()) {
			if (isNewLineChar((char)next)) {
				return;
			}
		}
	}


	protected static int readToken(char[] target, Reader reader) throws IOException {
		
		int input;
		int i;
		for (i = 0; i < target.length; i++) {
			input = reader.read();
	
			// end of file reached
			if (input == -1) {
				if (i == 0) {
					return -1;
				}
				
				break;
			} 
			
			// end of line
			else if (isNewLineChar((char)input)) {
				break;
			}
	
			// store char
			else {
				target[i] = (char) input;
			}
		}
		return i;
	}

	protected static boolean isNewLineChar(char c) {
		if (c == '\n' || c == '\r') {
			return true;
		} else {
			return false;
		}
	}
	
}


