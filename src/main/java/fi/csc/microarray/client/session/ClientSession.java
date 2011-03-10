package fi.csc.microarray.client.session;

import java.io.File;
import java.io.IOException;
import java.util.zip.ZipException;
import java.util.zip.ZipFile;

import org.apache.log4j.Logger;


public class ClientSession {

	public static final String SESSION_FILE_EXTENSION = "zip";
	public static final String SESSION_METADATA_FILENAME = "session.xml";
	
	private static final Logger logger = Logger.getLogger(ClientSession.class);
	
	
	public static boolean isValidSessionFile(File file) {
		if (file == null) {
			return false;
		}
		
		// does the file exist?
		if (!file.exists()) {
			return false;
		}

		// correct extension?
		if (!file.getName().endsWith("." + SESSION_FILE_EXTENSION)) {
			return false;
		}

		// is it a zip?
		ZipFile zipFile = null;
		try {
			try {
				 zipFile= new ZipFile(file);
			} catch (ZipException e) {
				return false;
			} catch (IOException e) {
				return false;
			}

			// does it contain the session metadata file
			if (zipFile.getEntry(SESSION_METADATA_FILENAME) == null) {
				return false;
			}
		
		} finally {
			if (zipFile != null) {
				try {
					zipFile.close();
				} catch (IOException e) {
					logger.warn("could not close zip file: " + file.getName());
				}
			}
		}

		return true;
	}
	
}
