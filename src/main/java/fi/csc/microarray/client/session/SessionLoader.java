package fi.csc.microarray.client.session;

import java.io.File;
import java.io.IOException;
import java.util.zip.ZipException;

import javax.xml.bind.JAXBException;
import javax.xml.parsers.ParserConfigurationException;

import org.w3c.dom.Document;
import org.xml.sax.SAXException;

import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.XmlUtil;

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
	
	
	public void loadSession() throws ZipException, IOException, JAXBException, SAXException, ParserConfigurationException {
		
		Document doc = XmlUtil.parseFile(sessionFile);
		
		String version = doc.getDocumentElement().getAttribute("format-version");
		
		if ("1".equals(version)) {
			// old format, use old loader
			SessionLoaderImpl1 impl = new SessionLoaderImpl1(sessionFile, dataManager, isDatalessSession);
			impl.loadSession();
			
		} else {
			// use new loader
			SessionLoaderImpl2 impl = new SessionLoaderImpl2(sessionFile, dataManager);
			impl.loadSession();
		}
	}

}
