package fi.csc.microarray.auth;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import javax.security.auth.callback.Callback;
import javax.security.auth.callback.CallbackHandler;
import javax.security.auth.callback.NameCallback;
import javax.security.auth.callback.PasswordCallback;
import javax.security.auth.callback.UnsupportedCallbackException;
import javax.security.auth.login.FailedLoginException;
import javax.security.auth.login.LoginContext;
import javax.security.auth.login.LoginException;

import org.apache.log4j.Logger;

import fi.csc.microarray.config.DirectoryLayout;

public class JaasAuthenticationProvider implements AuthenticationProvider {

	private static final String LOGIN_CONTEXT_NAME = "Chipster"; // login context name in JAAS configuration file
	private static final String CONFIG_FILE = "jaas.config";
	
	private static final Logger logger = Logger.getLogger(JaasAuthenticationProvider.class);
	
	public JaasAuthenticationProvider() throws IOException {
		initialize(true);
	}
	
	public JaasAuthenticationProvider(boolean legacyConfiguration) throws IOException {
		initialize(legacyConfiguration);
	}
	
	public boolean authenticate(String username, char[] password) {

		// get the login context
		LoginContext lc = null;
		try {
			lc = new LoginContext(LOGIN_CONTEXT_NAME, new SimpleCallbackHandler(username, password));
		} catch (LoginException le) {
			logger.error("Cannot create LoginContext. ", le);
			return false;
		} catch (SecurityException se) {
			logger.error("Cannot create LoginContext. ", se);
			return false;
		} 

		// authenticate
		try {
			// attempt authentication
			lc.login();
			// if we return with no exception, authentication succeeded
			logger.info("Authentication successful for " + username);
		} catch (FailedLoginException fle) {
			logger.info("Authentication failed for " + username);
			return false;
		} catch (LoginException le) {
			logger.error("Could not perform authentication for " + username, le);
			return false;
		}
		
		// authentication ok;
		return true;
	}

	
	private class SimpleCallbackHandler implements CallbackHandler {
		
		private String username;
		private char[] password;
		
		public SimpleCallbackHandler(String username, char[] password) {
			this.username = username;
			this.password = password;
		}

		public void handle(Callback[] callbacks) throws IOException,
		UnsupportedCallbackException {

			for (int i = 0; i < callbacks.length; i++) {
				if (callbacks[i] instanceof NameCallback) {
					NameCallback nc = (NameCallback)callbacks[i];
					nc.setName(this.username);

				} else if (callbacks[i] instanceof PasswordCallback) {
					PasswordCallback pc = (PasswordCallback)callbacks[i];
					pc.setPassword(this.password);
					
				} else {
					throw new UnsupportedCallbackException
					(callbacks[i], "Unrecognized Callback");
				}
			}
		}
	}
	
	
	/**
	 * If the jaas.config file does not exist in the working directory, create it using
	 * the jaas.config.default.
	 * 
	 * Set the java.security.auth.loing.config property as the path to the workdir/jaas.config.
	 * @param b 
	 * 
	 * @throws IOException
	 */
	private void initialize(boolean legacyConfiguration) throws IOException {
		
		File jaasConfigFile;
		if (legacyConfiguration) {
			jaasConfigFile = new File(DirectoryLayout.getInstance().getConfDir() + File.separator + CONFIG_FILE);
		} else {
			jaasConfigFile = new File("conf" + File.separator + CONFIG_FILE);
		}
		
		// if config file does not exist in the work dir, create it using the defaults
		if (!jaasConfigFile.exists()) {
			logger.info("JAAS config file " + jaasConfigFile + " does not exist, JAAS authentication provider failed.");
			throw new FileNotFoundException(jaasConfigFile + " not found");			
		}

		// set location of the config
		System.setProperty("java.security.auth.login.config", jaasConfigFile.getPath());
	}
}
