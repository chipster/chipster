package fi.csc.microarray.filebroker;

import it.sauronsoftware.cron4j.Scheduler;

import java.io.File;
import java.io.FilenameFilter;
import java.io.IOException;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Timestamp;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.Date;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.TimerTask;

import org.apache.commons.io.FileUtils;
import org.apache.log4j.Logger;
import org.joda.time.DateTime;
import org.joda.time.format.DateTimeFormatter;
import org.joda.time.format.ISODateTimeFormat;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.config.DirectoryLayout;

/**
 * Metadata server keeps track of files and sessions that are saved to long term storage space of the file broker.
 * It keeps track of metadata that is related to data files and session files. Runs on top of embedded Derby SQL 
 * database.
 * 
 * @author Aleksi Kallio
 *
 */
public class DerbyMetadataServer {
	
	private static final String DB_ROOT = "db-root";
	private static final String DB_NAME = "ChipsterFilebrokerMetadataDatabase";

	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(DerbyMetadataServer.class);
	
	private static final String METADATA_BACKUP_PREFIX="filebroker-metadata-db-backup-";

	public static final String DEFAULT_EXAMPLE_SESSION_OWNER = "example_session_owner";
	public static final String DEFAULT_EXAMPLE_SESSION_FOLDER = "Example sessions";

	
	private static String[][] SQL_CREATE_TABLES = new String[][] {
		{ 
			"sessions",
			"CREATE TABLE chipster.sessions (" +
					"uuid VARCHAR(200) PRIMARY KEY,  " +
					"name VARCHAR(200),  " +
				"username VARCHAR(200))" 
		},
		{
			"files",
			"CREATE TABLE chipster.files (" +
					"uuid VARCHAR(200) PRIMARY KEY,  " +
					"size BIGINT,  " +
					"created TIMESTAMP,  " +
					"last_accessed TIMESTAMP)"
		},
		{
			"belongs_to",
			"CREATE TABLE chipster.belongs_to (" + 
					"session_uuid VARCHAR(200)," +
					"file_uuid VARCHAR(200))"
		},
		{
			"special_users",
			"CREATE TABLE chipster.special_users (" + 
					"username VARCHAR(200) PRIMARY KEY," + 
					"show_as_folder VARCHAR(200))"
		}
	};

	
	private static String SQL_INSERT_SESSION  = "INSERT INTO chipster.sessions (name, username, uuid) VALUES (?, ?, ?)";
	private static String SQL_SELECT_SESSIONS_BY_USERNAME = "SELECT name, CAST(null AS VARCHAR(200)) as folder, uuid FROM CHIPSTER.SESSIONS WHERE username = ? UNION SELECT name, show_as_folder as folder, uuid FROM CHIPSTER.SESSIONS,  CHIPSTER.SPECIAL_USERS WHERE CHIPSTER.SESSIONS.username = CHIPSTER.SPECIAL_USERS.username";
	private static String SQL_SELECT_SESSIONS_BY_NAME_AND_USERNAME = "SELECT uuid FROM chipster.sessions WHERE name = ? AND username = ?";
	private static String SQL_DELETE_SESSION  = "DELETE FROM chipster.sessions WHERE uuid = ?";
	private static String SQL_SELECT_SESSIONS_BY_USERNAME_AND_UUID  = "SELECT uuid FROM chipster.sessions WHERE username = ? AND uuid = ?";
	private static String SQL_UPDATE_SESSION_NAME  = "UPDATE chipster.sessions SET name = ? WHERE uuid = ?";
	
	private static String SQL_INSERT_FILE  = "INSERT INTO chipster.files (uuid, size, created, last_accessed) VALUES (?, ?, ?, ?)";
	private static String SQL_UPDATE_FILE_ACCESSED  = "UPDATE chipster.files SET last_accessed = ? WHERE uuid = ?";
	private static String SQL_SELECT_FILE_BY_UUID = "SELECT * FROM chipster.files WHERE uuid = ?";
	private static String SQL_SELECT_FILES_TO_BE_ORPHANED  = "SELECT uuid from chipster.files WHERE uuid IN (SELECT file_uuid from chipster.belongs_to WHERE session_uuid = ?) AND uuid NOT IN (SELECT file_uuid from chipster.belongs_to WHERE NOT session_uuid = ?)";
	private static String SQL_DELETE_FILE  = "DELETE FROM chipster.files WHERE uuid = ?";
	
	private static String SQL_INSERT_BELONGS_TO  = "INSERT INTO chipster.belongs_to (session_uuid, file_uuid) VALUES (?, ?)";
	private static String SQL_DELETE_BELONGS_TO  = "DELETE FROM chipster.belongs_to WHERE session_uuid = ?";
	
	private static String SQL_INSERT_SPECIAL_USER  = "INSERT INTO chipster.special_users (username, show_as_folder) VALUES (?, ?)";
	// a sum of all distinct files referenced by the user's sessions
	// when multiple sessions contain the same files, their size is counted only once
	private static String SQL_LIST_STORAGE_USAGE_OF_USER = "SELECT SUM(size) as size FROM(SELECT DISTINCT chipster.files.uuid, chipster.files.size as size FROM chipster.sessions JOIN chipster.belongs_to ON chipster.sessions.uuid = chipster.belongs_to.session_uuid JOIN chipster.files ON chipster.files.uuid = chipster.belongs_to.file_uuid WHERE chipster.sessions.username = ?) as foo";
	// a sum of all distinct files referenced by the user's sessions
	// duplicate files are counted once for each user
	private static String SQL_LIST_STORAGE_USAGE_OF_USERS = "SELECT username, SUM(size) as size FROM (SELECT DISTINCT username, file_uuid, size FROM chipster.sessions JOIN chipster.belongs_to ON chipster.sessions.uuid = chipster.belongs_to.session_uuid JOIN chipster.files ON chipster.files.uuid = chipster.belongs_to.file_uuid) as foo GROUP BY username";
	// a simple sum of all files
	private static String SQL_LIST_STORAGE_USAGE_OF_SESSIONS = "SELECT chipster.sessions.username, chipster.sessions.name, chipster.sessions.uuid, SUM(chipster.files.size) AS size , MAX(chipster.files.last_accessed) AS date FROM chipster.sessions JOIN chipster.belongs_to ON chipster.sessions.uuid = chipster.belongs_to.session_uuid  JOIN chipster.files ON chipster.files.uuid = chipster.belongs_to.file_uuid WHERE chipster.sessions.username = ? GROUP BY chipster.sessions.uuid, chipster.sessions.name, chipster.sessions.username";
	// no join with sessions, so each file is counted only once (even when referenced by different users)
	private static String SQL_GET_TOTAL_DISK_USAGE = "SELECT SUM(chipster.files.size) AS size FROM chipster.files";
		
	private static String SQL_LIST_ALL_FILES = "SELECT * FROM chipster.files";

	// some statistics	
	private static String SQL_FILES_COUNT = "SELECT COUNT(*) FROM chipster.files";
	private static String SQL_SESSIONS_COUNT = "SELECT COUNT(*) FROM chipster.sessions";
	private static String SQL_MAPPINGS_COUNT = "SELECT COUNT(*) FROM chipster.belongs_to";
	private static String SQL_SPECIAL_USERS_COUNT = "SELECT COUNT(*) FROM chipster.special_users";
	// these result always 0 unless the database is corrupted 
	private static String SQL_ORPHAN_FILES = "SELECT COUNT(*) FROM chipster.files WHERE uuid NOT IN (SELECT session_uuid FROM chipster.belongs_to) AND uuid NOT IN (SELECT file_uuid FROM chipster.belongs_to)";
	private static String SQL_MISSING_FILES = "SELECT COUNT(*) FROM chipster.belongs_to WHERE file_UUID NOT IN (SELECT uuid FROM chipster.files)";
	private static String SQL_ORPHAN_SESSIONS = "SELECT COUNT(*) FROM chipster.sessions WHERE uuid NOT IN (SELECT session_uuid FROM chipster.belongs_to)";
	private static String SQL_MISSING_SESSIONS = "SELECT COUNT(*) FROM chipster.belongs_to WHERE session_UUID NOT IN (SELECT uuid FROM chipster.sessions)";

	private static String SQL_BACKUP = "CALL SYSCS_UTIL.SYSCS_BACKUP_DATABASE(?)";
	
	private Connection connection = null;

	/**
	 * Initialises the server. If underlying embedded Derby SQL database is not initialised, it is 
	 * initialised first.
	 * 
	 * @throws InstantiationException
	 * @throws IllegalAccessException
	 * @throws ClassNotFoundException
	 * @throws SQLException
	 * @throws IllegalConfigurationException 
	 * @throws IOException 
	 */
	public DerbyMetadataServer() throws InstantiationException, IllegalAccessException, ClassNotFoundException, SQLException, IOException, IllegalConfigurationException {
		
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
		
		// initialise connection
		System.setProperty("derby.system.home", DB_ROOT);
		Class.forName("org.apache.derby.jdbc.EmbeddedDriver").newInstance(); // allows multiple connections in one JVM, but not from multiple JVM's
		
		String strUrl = "jdbc:derby:" + DB_NAME + ";";
		
		File metadataBackupDir = DirectoryLayout.getInstance().getFilebrokerMetadataBackupDir();
		String restorePath = configuration.getString("filebroker", "metadata-restore-path");
		String fullRestorePath = metadataBackupDir + File.separator + restorePath + File.separator + DB_NAME;
		
		if (restorePath != null && !restorePath.isEmpty()) {
			logger.info("restoring metadata database from " + fullRestorePath);
			File dbDir = new File(DB_ROOT, DB_NAME);
			if (dbDir.exists()) {
				throw new IllegalConfigurationException("metadata restore isn't allowed, because the database " + dbDir + " exists already");
			}
			// One way to restore a backup
			// See http://db.apache.org/derby/docs/10.1/adminguide/tadminhubbkup44.html
			strUrl += "restoreFrom=" + fullRestorePath;
		} else {
			strUrl += "create=true";
		}
			
		
		connection = DriverManager.getConnection(strUrl);
		
		// initialise database, if needed
		initialise();

		logger.info("metadata database started");
		
		
		// initialise metadata database backup
		if (configuration.getBoolean("filebroker", "enable-metadata-backups")) {
		
			// get backup configuration
			int metadataBackupKeepCount = configuration.getInt("filebroker", "metadata-backup-keep-count");
			String backupTime = configuration.getString("filebroker", "metadata-backup-time").trim();
			
			// schedule backup tasks
			Scheduler scheduler = new Scheduler();
			scheduler.schedule(backupTime, new MetadataBackupTimerTask(metadataBackupDir, metadataBackupKeepCount));
			scheduler.start();
			
			logger.info("metadata backups enabled at: " + backupTime);
		} else {
			logger.info("metadata backups disabled");
		}	
	}
	
	private void initialise() throws SQLException {

		// create all missing tables
		int tableCount = 0;
		for (int i = 0; i < SQL_CREATE_TABLES.length; i++) {
			
			String table = SQL_CREATE_TABLES[i][0];
		
			ResultSet tables = connection.getMetaData().getTables(null, "CHIPSTER", table.toUpperCase(), new String[] { "TABLE" });
			if (!tables.next()) {
				
				// table does not exist, create it
				String createTable = SQL_CREATE_TABLES[i][1];
				PreparedStatement ps = connection.prepareStatement(createTable);
				ps.execute();
				tableCount++;
				
				// populate table, if needed
				if (table.equals("special_users")) {
					addSpecialUser(DEFAULT_EXAMPLE_SESSION_OWNER, DEFAULT_EXAMPLE_SESSION_FOLDER);
				}
			}			
		}
		
		// report what was done
		if (tableCount > 0) {
			logger.info("Created " + tableCount + " missing tables to database");
		}
	}
	
	/**
	 * List sessions that are available to all users, i.e. the example sessions.
	 * 
	 * @return
	 * @throws SQLException
	 */
	public List<DbSession> listPublicSessions() throws SQLException {
		return listSessions(null);
	}

	/**
	 * Lists all sessions that are available to the given user.
	 * 
	 * @param username
	 * @return
	 * @throws SQLException
	 */
	public List<DbSession> listSessions(String username) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_SELECT_SESSIONS_BY_USERNAME);
		ps.setString(1, username);
		ResultSet rs = ps.executeQuery();
		
		LinkedList<DbSession> sessions = new LinkedList<>();
		
		// go through files and add them, creating folders when needed
		HashSet<String> folders = new HashSet<>();
		while (rs.next()) {
			String name = rs.getString("name");
			if (rs.getString("folder") != null) {
				
				// need to make this file show inside a folder
				String folder = rs.getString("folder");
				
				// folder not yet seen, make entry for it first
				if (!folders.contains(folder)) {
					folders.add(folder);
					DbSession session = new DbSession("", folder + "/", username);
					sessions.add(session);
				}
				
				// prefix file name with folder
				name =  folder + "/" + name;
			}
			DbSession session = new DbSession(rs.getString("uuid"), name, username);
			sessions.add(session);
		}

		return sessions;
	}
	
	public DbFile fetchFile(String uuid) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_SELECT_FILE_BY_UUID);
		ps.setString(1, uuid);
		ResultSet rs = ps.executeQuery();
		
		if (rs.next()) {
			DbFile file = new DbFile(rs.getString(1), Long.parseLong(rs.getString(2)), rs.getString(3), rs.getString(4));
			return file;
		} else {
			return null;
		}
	}


	/**
	 * 'Touches' file.
	 * 
	 * @param uuid
	 * @throws SQLException
	 */
	public void markFileAccessed(String uuid) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_UPDATE_FILE_ACCESSED);
		ps.setTimestamp(1, new Timestamp(new Date().getTime()));
		ps.setString(2, uuid);
		ps.execute();
	}
	
	/**
	 * Adds metadata of a data file to the database.
	 * Doesn't throw an exception if the add fails because of the
	 * identical existing entry.
	 * 
	 * @param uuid unique identifier (filename) of the file
	 * @param size file size in bytes
	 * 
	 * @throws SQLException
	 */
	public void addFile(String uuid, long size) throws SQLException {
	
		try {
			addFileImpl(uuid, size);				
		} catch (SQLException e) {
			// don't care about the exception if the entry exists already
			try {			
				DbFile file = fetchFile(uuid);			
				if (file == null) {						
					throw e;
				} else {
					if (uuid.equals(file.getUuid()) && size == file.getSize()) {
						logger.debug("addFile failed, but the entry exist alreadỵ. Consider this as succesful");
						return;
					}
				}
			} catch (SQLException e2) {
				logger.debug("addFile failed and consequent fetchFile failed also. " +
						"The exception of the addFile is thrown and the exception of the fetchFile is logged here", e2);
				throw e;
			}
		}
	}

	private void addFileImpl(String uuid, long size) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_INSERT_FILE);
		ps.setString(1, uuid);
		ps.setLong(2, size);
		Timestamp now = new Timestamp(new Date().getTime());
		ps.setTimestamp(3, now);
		ps.setTimestamp(4, now);
		ps.execute();
	}

	/**
	 * Adds special username to the database. All sessions owned by the special
	 * users are visible to everyone. They can be used to create shared
	 * example sessions etc.
	 * 
	 * @param username special username to add
	 * @throws SQLException
	 */
	public void addSpecialUser(String username, String showAsFolder) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_INSERT_SPECIAL_USER);
		ps.setString(1, username);
		ps.setString(2, showAsFolder);
		ps.execute();
	}
	
	/**
	 * Links data file to a session (file).
	 * 
	 * @param fileUuid identifier of the data file
	 * @param sessionUuid identifier of the session
	 * 
	 * @throws SQLException
	 */
	public void linkFileToSession(String fileUuid, String sessionUuid) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_INSERT_BELONGS_TO);
		ps.setString(1, sessionUuid);
		ps.setString(2, fileUuid);
		ps.execute();
		
	}
	
	/**
	 * Adds session to the database.
	 * 
	 * @param username owner of the session
	 * @param name human readable name of the session
	 * @param uuid identifier of the session
	 * 
	 * @throws SQLException
	 */
	public void addSession(String username, String name, String uuid) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_INSERT_SESSION);
		ps.setString(1, name);
		ps.setString(2, username);
		ps.setString(3, uuid);
		ps.execute();
	}
	
	public void renameSession(String newName, String uuid) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_UPDATE_SESSION_NAME);
		ps.setString(1, newName);
		ps.setString(2, uuid);
		ps.execute();
	}
	
	public String fetchSession(String username, String name)  throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_SELECT_SESSIONS_BY_NAME_AND_USERNAME);
		ps.setString(1, name);
		ps.setString(2, username);
		ResultSet sessions = ps.executeQuery();
		if (sessions.next()) {
			return sessions.getString(1);
		} else {
			return null;
		}
	}

	/**
	 * Checks if user owns the session and can remove it.
	 * @param username 
	 * @param sessionUuid
	 * @return true if operation is allowed
	 * @throws SQLException 
	 */
	public boolean isUsernameAllowedToRemoveSession(String username, String sessionUuid) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_SELECT_SESSIONS_BY_USERNAME_AND_UUID);
		ps.setString(1, username);		
		ps.setString(2, sessionUuid);
		ResultSet rs = ps.executeQuery();

		return rs.next(); // return true if result set is not empty
	}
	
	
 	/**
	 * Removes session and dependent entries from the database. Dependent entries
	 * include file-session of the removed session and all linked files that are not
	 * linked to any other sessions.
	 * 
	 * @param uuid
	 * @throws SQLException
	 */
	public List<String> removeSession(String uuid) throws SQLException {

		// collect removed files so that they can be removed also physically
		LinkedList<String> removed = new LinkedList<String>();
		
		// find data files that will orphaned (must be done before removing belongs_to)
		PreparedStatement selectPs = connection.prepareStatement(SQL_SELECT_FILES_TO_BE_ORPHANED);
		selectPs.setString(1, uuid);
		selectPs.setString(2, uuid);
		ResultSet uuidRs = selectPs.executeQuery();
		LinkedList<String> orphanUuids = new LinkedList<String>();
		while (uuidRs.next()) {
			orphanUuids.add(uuidRs.getString(1));
		}

		// remove session entry from db 
		// ("entry point" is removed first, so if something fails, broken session entry is not left behind)
		PreparedStatement sessionPs = connection.prepareStatement(SQL_DELETE_SESSION);
		sessionPs.setString(1, uuid);
		sessionPs.execute();

		// remove belongs_to entry from db
		PreparedStatement belongsToPs = connection.prepareStatement(SQL_DELETE_BELONGS_TO);
		belongsToPs.setString(1, uuid);
		belongsToPs.execute();

		// remove session file entry from db and add to list of removed files
		PreparedStatement sessionFilePs = connection.prepareStatement(SQL_DELETE_FILE);
		sessionFilePs.setString(1, uuid);
		sessionFilePs.execute();
		removed.add(uuid);
		
		
		// remove orphaned data file entries from db and add to list of removed files
		for (String orphanUuid : orphanUuids) {
			PreparedStatement dataFilePs = connection.prepareStatement(SQL_DELETE_FILE);
			dataFilePs.setString(1, orphanUuid);
			dataFilePs.execute();
			removed.add(orphanUuid);
		}

		return removed;
	}
	
	/**
	 * Backup the database using the online backup procedure.
	 * 
	 * During the interval the backup is running, the database can be read, but writes to the database are blocked.
	 * 
	 * @param backupDir directory which will contain the db backup
	 * @throws SQLException
	 */
	public void backup(String backupDir) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_BACKUP);
		ps.setString(1, backupDir.replace(File.separator, "/"));
		ps.execute();
	}
	
	
	/**
	 * TimerTask for running metadata database backup.
	 * 
	 * @author hupponen
	 *
	 */
	private class MetadataBackupTimerTask extends TimerTask {

		private File baseBackupDir;
		private int backupKeepCount;
		
		/**
		 * Backup result will be backupDir/filebroker-metadata-db-backup-yyyy-MM-dd_mm:ss/ChipsterFilebrokerMetadataDatabase
		 * 
		 * @param backupDir base directory for backups, individual backups will be subdirectories 
		 */
		public MetadataBackupTimerTask(File backupDir, int backupKeepCount) {
			this.baseBackupDir = backupDir;
			this.backupKeepCount = backupKeepCount;
		}
		
		
		@Override
		public void run() {
			logger.info("backing up metadata database");
			long startTime = System.currentTimeMillis();
			DateFormat df = new SimpleDateFormat("yyyy-MM-dd_mm:ss");
			
			String fileName = baseBackupDir.getAbsolutePath() + File.separator + METADATA_BACKUP_PREFIX + df.format(new Date());
			try {
				backup(fileName);
			} catch (SQLException e) {
				logger.error("backing up metadata database failed", e);
			}
			logger.info("metadata backup took " + (System.currentTimeMillis() - startTime) + " ms");

			
			// remove old backups
			if (backupKeepCount >= 0) {

				// get backup dirs
				File[] backupDirs = baseBackupDir.listFiles(new FilenameFilter() {

					@Override
					public boolean accept(File dir, String name) {
						if (name.startsWith(METADATA_BACKUP_PREFIX)) {
							return true;
						}
						return false;
					}
				});

				// need to delete old?
				if (backupDirs.length > backupKeepCount) {
					long deleteStartTime = System.currentTimeMillis();
					
					// sort according to file name
					Arrays.sort(backupDirs, new Comparator<File>() {
						@Override
						public int compare(File o1, File o2) {
							return o1.getName().compareTo(o2.getName());
						}
						
					});
	
					// delete oldest until at keep limit
					for (int i = 0; backupDirs.length - i > backupKeepCount; i++) {
						logger.info("deleting old metadata backup: " + backupDirs[i].getName());
						try {
							FileUtils.deleteDirectory(backupDirs[i]);
						} catch (IOException e) {
							logger.error("could not delete old metadata backup directory: " + backupDirs[i]);
						}
					}
					logger.info("deleting old metadata backups took " + (System.currentTimeMillis() - deleteStartTime) + " ms");
				}
			}
		}
	}

	public Long getStorageusageOfUser(String username) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_LIST_STORAGE_USAGE_OF_USER);
		ps.setString(1, username);
		ResultSet rs = ps.executeQuery();

		Long size = null;
		if (rs.next()) {
			String sizeString = rs.getString("size");
			if (sizeString != null) {
				size = Long.parseLong(sizeString);
			} else {
				size = 0l; //user doesn't have any sessions
			}
		}
	
		return size;
	}

	@SuppressWarnings("unchecked")
	public List<String>[] getStorageusageOfUsers() throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_LIST_STORAGE_USAGE_OF_USERS);

		ResultSet rs = ps.executeQuery();
		LinkedList<String> usernames = new LinkedList<String>();
		LinkedList<String> sizes = new LinkedList<String>();

		while (rs.next()) {
			String username = rs.getString("username");
			String size = rs.getString("size");
			usernames.add(username);
			sizes.add(size);
		}

		return new List[] { usernames, sizes };
	}

	@SuppressWarnings("unchecked")
	public List<String>[] getStorageUsageOfSessions(String username) throws SQLException {
		
		PreparedStatement ps = connection.prepareStatement(SQL_LIST_STORAGE_USAGE_OF_SESSIONS);
		ps.setString(1, username);
		ResultSet rs = ps.executeQuery();
		
		LinkedList<String> usernames = new LinkedList<String>();
		LinkedList<String> sessions = new LinkedList<String>();
		LinkedList<String> sizes = new LinkedList<String>();
		LinkedList<String> dates = new LinkedList<String>();
		LinkedList<String> ids = new LinkedList<String>();

		
		DateTimeFormatter dateTimeFormatter = ISODateTimeFormat.dateTime();
		
		while (rs.next()) {
			String user = rs.getString("username");
			String session = rs.getString("name");
			String size = rs.getString("size");
			String id = rs.getString("uuid");
			DateTime date = new DateTime(rs.getTimestamp("date"));
			usernames.add(user);
			sessions.add(session);
			sizes.add(size);
			ids.add(id);
			dates.add(dateTimeFormatter.print(date));
		}

		return new List[] { usernames, sessions, sizes, dates, ids };
	}
	
	public String getStorageUsageTotal() throws SQLException {
		
		PreparedStatement ps = connection.prepareStatement(SQL_GET_TOTAL_DISK_USAGE);
		ResultSet rs = ps.executeQuery();
		
		rs.next();
		String size = rs.getString("size");
		
		if (size == null) {
			// when db is empty
			size = "0";
		}

		return size;
	}

	public List<DbFile> listAllFiles() throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_LIST_ALL_FILES);
		ResultSet rs = ps.executeQuery();
		
		List<DbFile> files = new ArrayList<DbFile>();
		
		while (rs.next()) {
			DbFile file = new DbFile(rs.getString(1), Long.parseLong(rs.getString(2)), rs.getString(3), rs.getString(4));
			files.add(file);
		}
		
		return files;
	}
	
	@SuppressWarnings("unchecked")
	public List<String>[] getStatistics() throws SQLException {
		
		String[] queries = { 
				SQL_FILES_COUNT, 
				SQL_SESSIONS_COUNT, 
				SQL_MAPPINGS_COUNT, 
				SQL_SPECIAL_USERS_COUNT, 
				SQL_ORPHAN_FILES, 
				SQL_MISSING_FILES, 
				SQL_ORPHAN_SESSIONS, 
				SQL_MISSING_SESSIONS 
		};
		
		String[] names = { 
				"rows in files table                ", 
				"rows in session table              ", 
				"mappings between files and sessions", 
				"number of special users            ", 
				"orphan files                       ", 
				"missing files                      ", 
				"orphan sessions                    ", 
				"missing sessions                   " };
	
		List<String> resultNames = new ArrayList<String>();
		List<String> resultvalues = new ArrayList<String>();
		
		for (int i = 0; i < queries.length && i < names.length; i++) {
		
			PreparedStatement ps = connection.prepareStatement(queries[i]);
			ResultSet rs = ps.executeQuery();

			rs.next();
			
			resultNames.add(names[i]);
			resultvalues.add(rs.getString(1));
		}
		
		return new List[] { resultNames, resultvalues };
	}
}
