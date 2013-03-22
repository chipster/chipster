package fi.csc.microarray.filebroker;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Timestamp;
import java.util.Date;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.Logger;

/**
 * Metadata server keeps track of files and sessions that are saved to long term storage space of the file broker.
 * It keeps track of metadata that is related to data files and session files. Runs on top of embedded Derby SQL 
 * database.
 * 
 * @author Aleksi Kallio
 *
 */
public class DerbyMetadataServer {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(DerbyMetadataServer.class);
	
	private static final String DEFAULT_EXAMPLE_SESSION_OWNER = "example_session_owner";

	private static String SCHEMA = "chipster";
	private static String SESSION_DBTABLE = "sessions";
	private static String FILE_DBTABLE = "files";
	private static String BELONGS_TO_DBTABLE = "belongs_to";
	private static String SPECIAL_USERS_DBTABLE = "special_users";
	
	private static String[][] SQL_CREATE_TABLES = new String[][] {
		{ 
			SESSION_DBTABLE,
			"CREATE TABLE " + SCHEMA + ". " + SESSION_DBTABLE + " (" +
					"uuid VARCHAR(200) PRIMARY KEY,  " +
					"name VARCHAR(200),  " +
				"username VARCHAR(200))" 
		},
		{
			FILE_DBTABLE,
			"CREATE TABLE " + SCHEMA + ". " + FILE_DBTABLE + " (" +
					"uuid VARCHAR(200) PRIMARY KEY,  " +
					"size BIGINT,  " +
					"created TIMESTAMP,  " +
					"last_accessed TIMESTAMP)"
		},
		{
			BELONGS_TO_DBTABLE,
			"CREATE TABLE " + SCHEMA + ". " + BELONGS_TO_DBTABLE + " (" + 
					"session_uuid VARCHAR(200) CONSTRAINT session_foreign_key REFERENCES " + SCHEMA + ". " + SESSION_DBTABLE + "," +
					"file_uuid VARCHAR(200) CONSTRAINT file_foreign_key REFERENCES " + SCHEMA + ". " + FILE_DBTABLE + ")"
		},
		{
			SPECIAL_USERS_DBTABLE,
			"CREATE TABLE " + SCHEMA + ". " + SPECIAL_USERS_DBTABLE + " (" + 
					"username VARCHAR(200) PRIMARY KEY)"
		}
	};

	private static String SQL_INSERT_SESSION  = "INSERT INTO " + SCHEMA + "." + SESSION_DBTABLE + " (name, username, uuid) VALUES (?, ?, ?)";
	private static String SQL_SELECT_SESSIONS_BY_USERNAME = "SELECT name, uuid FROM " + SCHEMA + "." + SESSION_DBTABLE + " WHERE username = ? OR username in (SELECT username FROM " + SCHEMA + "." + SPECIAL_USERS_DBTABLE + ")";
	private static String SQL_DELETE_SESSION  = "DELETE FROM " + SCHEMA + "." + SESSION_DBTABLE + " WHERE uuid = ?";
	
	private static String SQL_INSERT_FILE  = "INSERT INTO " + SCHEMA + "." + FILE_DBTABLE + " (uuid, size, created, last_accessed) VALUES (?, ?, ?, ?)";
	private static String SQL_UPDATE_FILE_ACCESSED  = "UPDATE " + SCHEMA + "." + FILE_DBTABLE + " SET last_accessed = ? WHERE uuid = ?";
	private static String SQL_SELECT_FILES_TO_BE_ORPHANED  = "SELECT uuid from " + SCHEMA + "." + FILE_DBTABLE + " WHERE uuid IN (SELECT file_uuid from " + SCHEMA + "." + BELONGS_TO_DBTABLE + " WHERE session_uuid = ?) AND uuid NOT IN (SELECT file_uuid from " + SCHEMA + "." + BELONGS_TO_DBTABLE + " WHERE NOT session_uuid = ?)";
	private static String SQL_DELETE_FILE  = "DELETE FROM " + SCHEMA + "." + FILE_DBTABLE + " WHERE uuid = ?";
	
	private static String SQL_INSERT_BELONGS_TO  = "INSERT INTO " + SCHEMA + "." + BELONGS_TO_DBTABLE + " (session_uuid, file_uuid) VALUES (?, ?)";
	private static String SQL_DELETE_BELONGS_TO  = "DELETE FROM " + SCHEMA + "." + FILE_DBTABLE + " WHERE uuid = ?";
	
	private static String SQL_INSERT_SPECIAL_USER  = "INSERT INTO " + SCHEMA + "." + SPECIAL_USERS_DBTABLE + " (username) VALUES (?)";
	private static String SQL_DELETE_SPECIAL_USER  = "DELETE FROM " + SCHEMA + "." + SPECIAL_USERS_DBTABLE + " WHERE username = ?";

	private Connection connection = null;

	/**
	 * Initialises the server. If underlying embedded Derby SQL database is not initialised, it is 
	 * initialised first.
	 * 
	 * @throws InstantiationException
	 * @throws IllegalAccessException
	 * @throws ClassNotFoundException
	 * @throws SQLException
	 */
	public DerbyMetadataServer() throws InstantiationException, IllegalAccessException, ClassNotFoundException, SQLException {
		
		// initialise connection
		System.setProperty("derby.system.home", "db-root");
		Class.forName("org.apache.derby.jdbc.EmbeddedDriver").newInstance(); // allows multiple connections in one JVM, but not from multiple JVM's
		String strUrl = "jdbc:derby:ChipsterFilebrokerMetadataDatabase;create=true";
		connection = DriverManager.getConnection(strUrl);
		
		// initialise database, if needed
		initialise();

		logger.info("Database started");
	}
	
	private void initialise() throws SQLException {

		// create all missing tables
		int tableCount = 0;
		for (int i = 0; i < SQL_CREATE_TABLES.length; i++) {
			
			String table = SQL_CREATE_TABLES[i][0];
		
			ResultSet tables = connection.getMetaData().getTables(null, SCHEMA.toUpperCase(), table.toUpperCase(), new String[] { "TABLE" });
			if (!tables.next()) {
				
				// table does not exist, create it
				String createTable = SQL_CREATE_TABLES[i][1];
				PreparedStatement ps = connection.prepareStatement(createTable);
				ps.execute();
				tableCount++;
				
				// populate table, if needed
				if (table.equals(SPECIAL_USERS_DBTABLE)) {
					addSpecialUser(DEFAULT_EXAMPLE_SESSION_OWNER);
				}
			}
			
		}
		
		// report what was done
		if (tableCount > 0) {
			logger.info("Created " + tableCount + " missing tables to database");
		}
	}

	/**
	 * Lists all sessions that are available to the given user.
	 * 
	 * @param username
	 * @return
	 * @throws SQLException
	 */
	@SuppressWarnings("unchecked")
	public List<String>[] listSessions(String username) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_SELECT_SESSIONS_BY_USERNAME);
		ps.setString(1, username);
		ResultSet rs = ps.executeQuery();
		LinkedList<String> names = new LinkedList<String>();
		LinkedList<String> uuids = new LinkedList<String>();
		
		while (rs.next()) {
			names.add(rs.getString("name"));
			uuids.add(rs.getString("uuid"));
		}

		return new List[] { names, uuids };
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
	 * 
	 * @param uuid unique identifier (filename) of the file
	 * @param size file size in bytes
	 * 
	 * @throws SQLException
	 */
	public void addFile(String uuid, long size) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_INSERT_FILE);
		ps.setString(1, uuid);
		ps.setLong(2, size);
		Timestamp now = new Timestamp(new Date().getTime());
		ps.setTimestamp(3, now);
		ps.setTimestamp(4, now);
		ps.execute();
	}
	
	/**
	 * Removes a data file from the database.
	 * 
	 * @param uuid unique identifier (filename) of the file
	 * @throws SQLException
	 */
	public void removeFile(String uuid) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_DELETE_FILE);
		ps.setString(1, uuid);
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
	public void addSpecialUser(String username) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_INSERT_SPECIAL_USER);
		ps.setString(1, username);
		ps.execute();
	}
	
	/**
	 * Removes a special user from the database.
	 * 
	 * @param uuid special username to remove
	 * @throws SQLException
	 */
	public void removeSpecialUser(String username) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_DELETE_SPECIAL_USER);
		ps.setString(1, username);
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

	/**
	 * Removes session and dependent entries from the database. Dependent entries
	 * include file-session of the removed session and all linked files that are not
	 * linked to any other sessions.
	 * 
	 * @param uuid
	 * @throws SQLException
	 */
	public void removeSession(String uuid) throws SQLException {

		// unwind dependency structure: sessions, belongs_to, files
		// start by removing session
		PreparedStatement sessionPs = connection.prepareStatement(SQL_DELETE_SESSION);
		sessionPs.setString(1, uuid);
		sessionPs.execute();
		
		// find files that will be orphaned by removing this session
		PreparedStatement selectPs = connection.prepareStatement(SQL_SELECT_FILES_TO_BE_ORPHANED);
		selectPs.setString(1, uuid);
		ResultSet uuidRs = selectPs.executeQuery();
		LinkedList<String> orphanUuids = new LinkedList<String>();
		while (uuidRs.next()) {
			orphanUuids.add(uuidRs.getString(1));
		}
		
		// remove links
		PreparedStatement belongsToPs = connection.prepareStatement(SQL_DELETE_BELONGS_TO);
		belongsToPs.setString(1, uuid);
		belongsToPs.execute();
		
		// remove linked orphaned files
		for (String orphanUuid : orphanUuids) {
			PreparedStatement filePs = connection.prepareStatement(SQL_DELETE_FILE);
			filePs.setString(1, orphanUuid);
			filePs.execute();
		}
	}
	
	
}
