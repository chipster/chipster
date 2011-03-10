package fi.csc.microarray.client.session;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Reader;
import java.util.LinkedHashMap;
import java.util.zip.ZipFile;

import org.apache.log4j.Logger;
import org.w3c.dom.Document;

import fi.csc.microarray.databeans.DataFolder;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.XmlUtil;

public class SessionLoader {
	
	private File sessionFile;
	private Document sessionDoc;
	
	private LinkedHashMap<String, DataFolder> folders = new LinkedHashMap<String, DataFolder>();

	private static final Logger logger = Logger.getLogger(SessionLoader.class);
	
	
	public SessionLoader(File sessionFile) throws MicroarrayException {
		if (!ClientSession.isValidSessionFile(sessionFile)) {
			throw new MicroarrayException("Not a valid session file.");
		}
		this.sessionFile = sessionFile;
	}
	
	public void loadSession() {
		ZipFile zipFile = null;
		try {
			// get the session.xml zip entry
			zipFile = new ZipFile(sessionFile);
			Reader metadataReader = new BufferedReader(new InputStreamReader(zipFile.getInputStream(zipFile.getEntry(ClientSession.SESSION_METADATA_FILENAME))));

			// create the dom for session.xml
			this.sessionDoc = XmlUtil.parseReader(metadataReader);
			this.sessionElement = sessionDoc.getDocumentElement();
			
			
			parseFolders();
			
		} 
		// FIXME
		catch (Exception e) {
			logger.error(e);
		}

		// try to close all input streams from the zip file
		finally {
			if (zipFile != null) {
				try {
					zipFile.close();
				} catch (IOException e) {
					logger.warn("could not close zip file");
				}
			}
		}
	}

	private void parseFolders() {
		XmlUtil.getChildElements(sessionDoc., "folder");
		
	}
}
