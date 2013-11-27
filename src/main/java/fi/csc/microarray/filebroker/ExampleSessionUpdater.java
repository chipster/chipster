package fi.csc.microarray.filebroker;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.sql.Timestamp;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import de.schlichtherle.truezip.file.TFile;

/**
 * Synchronize example sessions in zip format and stored as server sessions.
 * 
 * Storing example session as server sessions enables many handy features, like deduplication 
 * of the derived datasets and easy modification of these session in the client. On the contrary, example
 * sessions stored as zip files are easy to handle in update scripts. This class tries to take 
 * advantages of both storage methods by keeping example session stored both as server sessions 
 * and zip files.
 * 
 * It is assumed that example sessions of the running instance are modified in the client and example
 * sessions of the stopped instance are modified in zip files. When the FileServer starts, the 
 * server sessions are updated to the same state with zip files. When the server sessions are modified
 * on the fly, the same operations are also applied to zip sessions.
 * 
 * @author klemela
 */
public class ExampleSessionUpdater extends FileServerListener {
	
	/* tar is fast and the content is compressed already (for example saving 40MB takes:
	 * - 300ms for tar
	 * - 700ms for zip
	 * - 2000ms for tar.gz 
	 */
	private static final String EXAMPLE_SESSION_ARCHIVE = "all-example-sessions.tar";

	private int ZIP_TIMESTAMP_PENALTY = 2000; //ms
	
	private static Logger logger = Logger.getLogger(ExampleSessionUpdater.class);

	private File exampleSessionDir;
	public DerbyMetadataServer metadataServer;
	private ServerSessionImportExportTool importExportTool;

	private FileServer fileServer;

	public ExampleSessionUpdater(FileServer fileServer, DerbyMetadataServer metadataServer, File exampleSessionDir) throws Exception {
		this.fileServer = fileServer;
		this.metadataServer = metadataServer;
		this.exampleSessionDir = exampleSessionDir;
		importExportTool = new ServerSessionImportExportTool(fileServer);
		
		if (!this.exampleSessionDir.exists()) {
			if (!this.exampleSessionDir.mkdir()) {
				throw new FileBrokerException("Creating example session directory failed: " + this.exampleSessionDir);
			}
		}
				
		fileServer.addListener(this);
	}

	/**
	 * Copy status of the zip sessions to server sessions:
	 * - new sessions are imported
	 * - existing sessions are updated if the zip is newer than the server session
	 * - server sessions without corresponding zip session are removed
	 * 
	 * @throws JMSException
	 * @throws SQLException
	 * @throws Exception
	 */
	public void importExampleSessions() throws JMSException, SQLException, Exception {	
		
		//list public server sessions
		//map keys are file name part of the server session name (called basename), while the DbSession.getName() contains also the path 
		List<DbSession> dbSessionList = metadataServer.listPublicSessions();    		
		Map<String, DbSession> dbSessions = new HashMap<>();    		    		
		for (DbSession session : dbSessionList) {
			String basename = session.getBasename();
			//filter out directories
			if (basename != null) {
				dbSessions.put(basename, session);
			}
		}
		//list zip file sessions
		//map keys are zip file names without their file extension (called basename)
		logger.debug("searching example sessions from " + exampleSessionDir);	
		HashMap<String, File> zipSessions = new HashMap<>();
		for (File file : exampleSessionDir.listFiles()) {
			if (isZipSession(file)) {				
				String basename = importExportTool.filenameToBasename(file.getName());
				zipSessions.put(basename, file);
			}
		}
		
		//if zip was removed, remove also the server session
		Iterator<String> dbSessionIter = dbSessions.keySet().iterator();
		while(dbSessionIter.hasNext()) {
			String dbSessionBasename = dbSessionIter.next();
			if (!zipSessions.containsKey(dbSessionBasename)) {
				logger.info("found a server session  '" + dbSessionBasename + "', but no zip session with that name. Going to remove the server session");
				fileServer.removeSession(dbSessions.get(dbSessionBasename).getUuid());
				dbSessionIter.remove();
			}
		}				
		
		//if a zip was added or updated, then import it						
		for (String zipSessionBasename : zipSessions.keySet()) {
			
			File zipSessionFile = zipSessions.get(zipSessionBasename);
			DbSession dbSession = dbSessions.get(zipSessionBasename);
			boolean store = false;
			
			if (dbSession == null) {
				logger.debug("found a zip session '" + zipSessionBasename + "', but there is no server session with that name");
				store = true;
			} else {
				//there is a server session with same name, check timestamps		
				DbFile dbSessionFile = metadataServer.fetchFile(dbSession.getUuid());
				Timestamp dbSessionCreated = Timestamp.valueOf(dbSessionFile.getCreated());
				//logger.debug("server session created " + dbSessionCreated);

				/* When a session is saved, its files are first upload and only after that the session is inserted into database.
				 * The session file gets its timestamp when the session file is uploaded, but we can export the session to a zip
				 * file only afterwards. 
				 * 
				 * As a consequence, it looks like the zip file is newer and needs to again imported, although there is no need for that.
				 * We have to make zip file to look little bit older than it really is in this comparison to avoid this extra import.
				 */
				Timestamp zipSessionModified = new Timestamp(zipSessionFile.lastModified() - ZIP_TIMESTAMP_PENALTY);
				//logger.debug("zip session modified: " + zipSessionModified);								
				
				store = zipSessionModified.after(dbSessionCreated);
			}

			if (store) {				
				//import the zip session
				logger.info("storing example session '" + zipSessionBasename + "'");
				importExportTool.importSession(zipSessionFile);
				logger.debug("example session stored");				
			} else {
				logger.debug("example session '" + zipSessionBasename + "' is up-to-date");
			}
		}
	}

	/**
	 * Export server session to zip session. Possible existing files are overwritten.
	 * 
	 * @param sessionUuid
	 * @param basename
	 */
	public void exportServerSession(String sessionUuid, String basename) {
		try {							
			File zipFile = basenameToFile(basename);
			removeIfExists(zipFile);
			
			importExportTool.exportSession(sessionUuid, zipFile);
			
			updateExampleSessionArchive();
			
		} catch (Exception e) {
			logger.error("session export failed", e);
		}
	}
	
	/**
	 * Store all session files to a tarball which is easy to copy when the example sessions of this 
	 * installation are needed elsewhere.
	 * 
	 * @throws IOException
	 */
	private void updateExampleSessionArchive() throws IOException {
		logger.debug("updating the archive of all example sessions");

		File archive = new File(exampleSessionDir, EXAMPLE_SESSION_ARCHIVE);
		archive.delete();
		
		for (File file : exampleSessionDir.listFiles()) {
			if (isZipSession(file)) {			
				TFile inputSession = new TFile(file);
				TFile outputSession = new TFile(archive, file.getName());
				inputSession.cp_rp(outputSession);
			}
		}
		TFile.umount();
	}

	private boolean isZipSession(File file) {
		return file.getName().endsWith(".zip") && !EXAMPLE_SESSION_ARCHIVE.equals(file.getName());
	}

	/**
	 * Get a File object for the example session basename (session name without directories and file extension)
	 * 
	 * @param basename
	 * @return
	 */
	private File basenameToFile(String basename) {
		return new File(this.exampleSessionDir, importExportTool.basenameToFileName(basename));
	}
	
	/**
	 * Iterate through the example sessions to find a name for the session uuid.
	 * 
	 * @param uuid
	 * @return
	 * @throws SQLException
	 */
	private String uuidToBasename(String uuid) throws SQLException {
		for (DbSession session : metadataServer.listPublicSessions()) {
			if (uuid.equals(session.getUuid())) {
				return session.getBasename();
			}
		}
		return null;
	}

	/**
	 * Find a session name for the session uuid and remove a zip session of that name.
	 * 
	 * @param uuid
	 */
	public void removeZipSession(String uuid) {
		try {
			String basename = uuidToBasename(uuid);			
			File zipFile = basenameToFile(basename);
			removeIfExists(zipFile);
			
			updateExampleSessionArchive();
		} catch (SQLException | IOException e) {
			logger.error("removing example session failed ", e);
		}
	}
	
	private void removeIfExists(File file) {
		if (file.exists()) {
			file.delete();
		}
	}

	/* 
	 * When something is done for the server sessions, do the some for the zip sessions:
	 * - if a new server session is stored, export it to zip session
	 * - if a server session is removed, remove also the zip session
	 * 
	 * Client will show an error message if this takes longer than a few seconds, because the
	 * session saving is just a database write for normal users  and client has very short 
	 * timeout for it.
	 * 
	 * (non-Javadoc)
	 * @see fi.csc.microarray.filebroker.FileServerListener#listen(fi.csc.microarray.filebroker.FileServerListener.Event)
	 */
	@Override
	public void listen(Event e) {
		
		//don't care about events originated in file broker (for example when we import a zip file, don't export it again here)
		if (!(e.getEndpoint() instanceof DirectFileBrokerEndpoint)) {
			
			if (e instanceof BeforeStoreSession) {
				BeforeStoreSession event = (BeforeStoreSession) e;

				if (DerbyMetadataServer.DEFAULT_EXAMPLE_SESSION_OWNER.equals(event.getUsername())) {
					logger.info("example session " + event.getSessionName() + " is being saved, exporting it also as a zip file");
					exportServerSession(event.getUuid(), event.getSessionName());
				}			
			}

			if (e instanceof BeforeRemoveSession) {
				BeforeRemoveSession event = (BeforeRemoveSession) e;

				if (DerbyMetadataServer.DEFAULT_EXAMPLE_SESSION_OWNER.equals(event.getUsername())) {
					logger.info("example session is being removed, removing also the zip file");
					removeZipSession(event.getUuid());
				}			
			}
		}
	}
}
