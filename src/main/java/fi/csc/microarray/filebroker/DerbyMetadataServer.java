package fi.csc.microarray.filebroker;

import org.apache.log4j.Logger;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.LinkedList;
import java.util.List;

public class DerbyMetadataServer {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(DerbyMetadataServer.class);

	private Connection connection = null;

	private static String SCHEMA = "chipster";
	private static String SESSION_DBTABLE = "sessions";
	private static String FILE_DBTABLE = "files";
	private static String BELONGS_TO_DBTABLE = "belongs_to";
	private static String SPECIAL_USERS_DBTABLE = "special_users";
	
	private static String[][] SQL_CREATE_TABLES = new String[][] {
		{ 
			"sessions",
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

	private static String SQL_SELECT_SESSIONS_BY_USERNAME = "SELECT name, uuid FROM " + SCHEMA + "." + SESSION_DBTABLE + " WHERE username = ? OR username in (SELECT username FROM " + SCHEMA + "." + SPECIAL_USERS_DBTABLE + ")";
	private static String SQL_INSERT_SESSION  = "INSERT INTO " + SCHEMA + "." + SESSION_DBTABLE + " (name, username, uuid) VALUES (?, ?, ?)";
	private static String SQL_DELETE_SESSION  = "DELETE FROM " + SCHEMA + "." + SESSION_DBTABLE + " WHERE uuid = ?";
	
	
	public DerbyMetadataServer() throws InstantiationException, IllegalAccessException, ClassNotFoundException, SQLException {
		
		// initialise connection
		System.setProperty("derby.system.home", "db-root");
		Class.forName("org.apache.derby.jdbc.EmbeddedDriver").newInstance(); // allows multiple connections in one JVM, but not from multiple JVM's
		String strUrl = "jdbc:derby:MetadataDatabase;create=true";
		connection = DriverManager.getConnection(strUrl);
		
		// initialise database, if needed
		initialise();

		logger.info("Database started");
	}
	
	private void initialise() throws SQLException {
		
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
			}
		}
		
		if (tableCount > 0) {
			logger.info("Created " + tableCount + " missing tables to database");
		}
	}

	@SuppressWarnings("unchecked")
	public List<String>[] listSessionsInDatabase(String username) throws SQLException {
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

	public void addSessionToDatabase(String username, String name, String uuid) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_INSERT_SESSION);
		ps.setString(1, name);
		ps.setString(2, username);
		ps.setString(3, uuid);
		ps.execute();
	}

	public void removeSessionFromDatabase(String uuid) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_DELETE_SESSION);
		ps.setString(1, uuid);
		ps.execute();
	}
	
	
}
