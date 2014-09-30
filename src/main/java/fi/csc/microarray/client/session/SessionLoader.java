package fi.csc.microarray.client.session;

import java.io.File;
import java.io.InputStreamReader;
import java.util.List;

import org.w3c.dom.Document;

import de.schlichtherle.truezip.zip.ZipFile;
import fi.csc.microarray.client.operation.OperationRecord;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.XmlUtil;

public class SessionLoader {
	
	private DataManager dataManager;
	private File sessionFile;
	private String sessionId;
	private boolean isDatalessSession;
	private Integer xOffset;
	private String sessionNotes;
	private List<OperationRecord> unfinishedJobs;

	public SessionLoader(File sessionFile, boolean isDatalessSession, DataManager dataManager) throws MicroarrayException {
		this.sessionFile = sessionFile;
		this.sessionId = null;
		this.dataManager = dataManager; 
		this.isDatalessSession = isDatalessSession;
	}

	public SessionLoader(String sessionId, DataManager dataManager) throws MicroarrayException {
		this.sessionFile = null;
		this.sessionId = sessionId;
		this.dataManager = dataManager; 
		this.isDatalessSession = true;
	}

	
	/**
	 * @return a list of OperationRecords for tasks that were running when the
	 *         session was saved
	 *         
	 * @throws Exception
	 */
	public void loadSession() throws Exception {
		
		ZipFile zipFile = null;
		InputStreamReader metadataReader = null;
		String version = Integer.toString(UserSession.SESSION_VERSION);
		try {
			// get the session.xml zip entry (only if a file, remote sessions are always latest version)
			if (sessionFile != null) {
				zipFile = new ZipFile(sessionFile);
				metadataReader = new InputStreamReader(zipFile.getInputStream(zipFile.getEntry(UserSession.SESSION_DATA_FILENAME)));
				Document doc = XmlUtil.parseReader(metadataReader);
				version = doc.getDocumentElement().getAttribute("format-version");
			}

		} finally {
			IOUtils.closeIfPossible(metadataReader);
			IOUtils.closeIfPossible(zipFile);
		}
		
		if (Integer.toString(UserSession.PREVIOUS_SESSION_VERSION).equals(version)) {
			// old format, use old loader
			SessionLoaderImpl1 impl = new SessionLoaderImpl1(sessionFile, dataManager, isDatalessSession);
			impl.loadSession();
			return;
			
		} else {
			// use new loader
			SessionLoaderImpl2 impl;
			if (sessionFile != null) {
				impl = new SessionLoaderImpl2(sessionFile, dataManager, isDatalessSession);
			} else {
				impl = new SessionLoaderImpl2(sessionId, dataManager, isDatalessSession);
			}
			impl.setXOffset(xOffset);
			impl.loadSession();
			sessionNotes = impl.getSessionNotes();
			unfinishedJobs = impl.getUnfinishedOperations();
		}
	}

	public void setXOffset(Integer xOffset) {
		this.xOffset = xOffset;
	}

	public String getSessionNotes() {
		return this.sessionNotes;
	}

	public List<OperationRecord> getUnifinishedJobs() {
		return unfinishedJobs;
	}
}
