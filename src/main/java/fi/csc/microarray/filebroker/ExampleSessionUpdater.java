package fi.csc.microarray.filebroker;

import java.io.File;
import java.sql.SQLException;
import java.sql.Timestamp;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.DirectMessagingEndpoint;

/**
 * Import example sessions in zip format and store them as server sessions.
 * 
 * Storing example session as server sessions enables many handy features, like
 * deduplication of the datasets and easy modification of these session in the
 * client. On the contrary, example sessions stored as zip files are easy to
 * handle in server update scripts.
 *
 * When the FileServer starts, the server sessions are updated to the
 * same state with zip files.
 * 
 * @author klemela
 */
public class ExampleSessionUpdater {
	
	private static Logger logger = Logger.getLogger(ExampleSessionUpdater.class);

	private File exampleSessionDir;
	public DerbyMetadataServer metadataServer;
	private ServerSessionImportExportTool importExportTool;

	private FileServer fileServer;

	public ExampleSessionUpdater(FileServer fileServer, DerbyMetadataServer metadataServer, File exampleSessionDir) throws Exception {
		this.fileServer = fileServer;
		this.metadataServer = metadataServer;
		this.exampleSessionDir = exampleSessionDir;
		
		//set up JMS replacement
		DirectMessagingEndpoint directEndpoint = new DirectMessagingEndpoint(DerbyMetadataServer.DEFAULT_EXAMPLE_SESSION_OWNER);
		fileServer.addEndpoint(directEndpoint);
		
		importExportTool = new ServerSessionImportExportTool(directEndpoint);
		
		if (!this.exampleSessionDir.exists()) {
			logger.info("example session directory " + this.exampleSessionDir + " doesn't exist");
			logger.info("example session import is disabled");
		}
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
		
		if (this.exampleSessionDir.exists()) {

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
					fileServer.removeSession(dbSessions.get(dbSessionBasename).getDataId());
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
					Timestamp dbSessionCreated = getTimestamp(dbSession.getDataId());				
					Timestamp zipSessionModified = getTimestamp(zipSessionFile);				

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
	}

	private Timestamp getTimestamp(File file) {
		return new Timestamp(file.lastModified());
	}

	private Timestamp getTimestamp(String dataId) throws SQLException {
		DbFile dbSessionFile = metadataServer.fetchFile(dataId);
		return Timestamp.valueOf(dbSessionFile.getCreated());
	}
		
	private boolean isZipSession(File file) {
		return file.getName().endsWith(".zip");
	}		
}
