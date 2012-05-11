package fi.csc.microarray.client.session;

import java.io.File;
import java.io.IOException;
import java.util.zip.ZipException;

import javax.xml.bind.JAXBException;

import org.xml.sax.SAXException;

import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.exception.MicroarrayException;

public class SessionLoader {
	
	private DataManager dataManager;
	private File sessionFile;
	private boolean isDatalessSession;

	public SessionLoader(File sessionFile, boolean isDatalessSession, DataManager dataManager) throws MicroarrayException {
		if (!UserSession.isValidSessionFile(sessionFile)) {
			throw new MicroarrayException("Not a valid session file.");
		}
		this.sessionFile = sessionFile;
		this.dataManager = dataManager; 
		this.isDatalessSession = isDatalessSession;
	}
	
	
	public void loadSession() throws ZipException, IOException, JAXBException, SAXException {

		SessionLoaderImpl1 impl = new SessionLoaderImpl1(sessionFile, dataManager, isDatalessSession);
		impl.loadSession();
	}

}
