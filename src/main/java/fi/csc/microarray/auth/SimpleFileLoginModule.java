package fi.csc.microarray.auth;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Map;

import javax.security.auth.Subject;
import javax.security.auth.callback.CallbackHandler;

import org.apache.log4j.Logger;
import org.mortbay.util.IO;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.MicroarrayConfiguration;
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.LookaheadStringReader;

/**
 * Login module for Chipster type user lists. They have format
 * username:password:expiration, where expiration is optional and there can by
 * anything after the actual content of the line. Also comments (#) can be used,
 * but they must start at first character of the line.
 * 
 * @author Taavi Hupponen, Aleksi Kallio
 * 
 */
public class SimpleFileLoginModule extends LoginModuleBase {

	public static final DateFormat EXPIRATION_DATE_FORMAT = new SimpleDateFormat("yyyy-MM-dd");
	public static final String DELIMETER_CHARACTER = ":";
	public static final String COMMENT_CHARACTER = "#";
	
	private static final Logger logger = Logger.getLogger(SimpleFileLoginModule.class);
	private static final String DEFAULT_PASSWD_FILE = "/users.default";

	// configurable options
	protected File passwdFile;

	public void initialize(Subject subject, CallbackHandler callbackHandler, Map<String, ?> sharedState, Map<String, ?> options) {
		super.initialize(subject, callbackHandler, sharedState, options);

		// passwd file

		// if file option is not defined or is empty, create the template file
		String passwdFileName = (String) options.get("passwdFile");
		if (passwdFileName == null || passwdFileName.equals("")) {
			createDefaultPasswdFile();
			return;
		}

		this.passwdFile = new File(passwdFileName);

		// absolute path
		if (passwdFile.isAbsolute()) {
			if (!passwdFile.exists()) {
				logger.warn("Passwd file " + passwdFile.getPath() + " not found.");
			}
		}

		// relative path, search the file from work dir, if not found create the
		// template
		else {
			passwdFile = new File(MicroarrayConfiguration.getWorkDir(), passwdFile.getName());
			if (!passwdFile.exists()) {
				logger.warn("Passwd file " + passwdFile.getPath() + " not found.");
				createDefaultPasswdFile();
			}
		}
	}

	protected boolean authenticate(String username, char[] password) throws IOException {

		logger.debug(this.getClass().getName() + " authenticating " + username);

		// passwd file
		BufferedReader reader = null;

		try {
			reader = new BufferedReader(new FileReader(this.passwdFile));
			// loop the lines of the password file
			for (String line = reader.readLine(); line != null; line = reader.readLine()) {

				// FIXME line should be cleaned (contains password)
				
				// check for empty line
				if (line.trim().length() == 0) {
					continue;
				}

				LookaheadStringReader tokens = new LookaheadStringReader(line);

				// check for comment line
				if (tokens.lookahead().equals(COMMENT_CHARACTER)) {
					continue;
				}

				// username
				String readUsername = tokens.readTo(DELIMETER_CHARACTER);
				if (!readUsername.equals(username)) {
					// did not match
					continue;
				}

				// delimiter (:)
				tokens.read();

				// password
				StringBuffer readPassword = tokens.readToSB(DELIMETER_CHARACTER);
				boolean match = true; 
				for (int i = 0; i < readPassword.length(); i++) {
					// check that we match (length checking is a bit redundant, but who cares...)
					if (readPassword.length() != password.length || readPassword.charAt(i) != password[i]) {
						match = false;
					}
					// clean password as we go
					readPassword.setCharAt(i, (char)0); 
				}

				if (!match) {
					// did not match
					continue;
				}

				// delimiter (:), if any
				if (!tokens.isAtEnd() && tokens.lookahead().equals(DELIMETER_CHARACTER)) {
					tokens.read();

					String expiration = tokens.readTo(DELIMETER_CHARACTER);
					
					if (expiration.trim().length() > 0) {
						// check only if data is not empty
						try {
							Date expirationDate = EXPIRATION_DATE_FORMAT.parse(expiration);
							if (new Date().after(expirationDate)) {
								match = false; // authentication successful, but account has expired
							}
						} catch (ParseException e) {
							logger.error("when authenticating " + username +" failed to parse exp. date: " + expiration);
							match = false;
						}
					}
				}

				if (match) {
					return true;
				}
			}

		} catch (Exception e) {
			e.printStackTrace();
			
		} finally {
			IOUtils.closeIfPossible(reader);
		}

		// matching line was not found
		return false;
	}

	private void createDefaultPasswdFile() {
		try {

			File defaultFile = new File(MicroarrayConfiguration.getWorkDir(), "users");
			if (defaultFile.exists()) {
				return;
			}

			logger.info("Creating new default passwd file " + defaultFile.getPath());
			InputStream defaults = Configuration.class.getResourceAsStream(DEFAULT_PASSWD_FILE);
			OutputStream out = new FileOutputStream(defaultFile);
			try {
				IO.copy(defaults, out);
			} finally {
				if (defaults != null) {
					defaults.close();
				}
				if (out != null) {
					out.close();
				}
			}
		} catch (IOException ioe) {
			logger.error("Could not create default users file.", ioe);
		}
	}

}
