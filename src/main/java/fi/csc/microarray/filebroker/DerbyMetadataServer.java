package fi.csc.microarray.filebroker;

import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.LinkedList;
import java.util.List;

public class DerbyMetadataServer {

	private Connection connection = null;

	private static String SESSION_DBTABLE = "sessions";
	private static String FILE_DBTABLE = "files";
	private static String BELONGS_TO_DBTABLE = "belongs_to";
	private static String[] DBTABLES = new String[] {
		SESSION_DBTABLE,
		FILE_DBTABLE,
		BELONGS_TO_DBTABLE
	};
	
	private static String[] SQL_CREATE_TABLES = new String[] {
		"CREATE TABLE " + SESSION_DBTABLE + " (" +
			"uuid VARCHAR(200) PRIMARY KEY,  " +
			"name VARCHAR(200),  " +
			"username VARCHAR(200))",
		"CREATE TABLE " + FILE_DBTABLE + " (" +
			"uuid VARCHAR(200) PRIMARY KEY,  " +
			"size BIGINT,  " +
			"created TIMESTAMP,  " +
			"last_accessed TIMESTAMP)",
		"CREATE TABLE " + BELONGS_TO_DBTABLE + " (" + 
			"session_uuid VARCHAR(200) CONSTRAINT session_foreign_key REFERENCES " + SESSION_DBTABLE + "," +
			"file_uuid VARCHAR(200) CONSTRAINT file_foreign_key REFERENCES " + FILE_DBTABLE + ")"
	};

	private static String SQL_SELECT_SESSIONS_BY_USERNAME = "SELECT * FROM " + SESSION_DBTABLE + " WHERE username = ?";
	private static String SQL_INSERT_SESSION  = "INSERT INTO " + SESSION_DBTABLE + " (name, username, uuid) VALUES (?, ?, ?)";
	private static String SQL_DELETE_SESSION  = "DELETE FROM " + SESSION_DBTABLE + " WHERE uuid = ?";
	
	
	public DerbyMetadataServer() throws InstantiationException, IllegalAccessException, ClassNotFoundException, SQLException {
		
		// initialise database
		System.setProperty("derby.system.home", "db-root");
		Class.forName("org.apache.derby.jdbc.EmbeddedDriver").newInstance(); // allows multiple connections in one JVM, but not from multiple JVM's
		String strUrl = "jdbc:derby:TestDatabase;create=true";
		connection = DriverManager.getConnection(strUrl);
		
		// initialise tables, if empty
		if (!isInitialised()) {
			initialise();
		}

		System.out.println("Database started");
	}
	
	private boolean isInitialised() throws SQLException {
		for (String table : DBTABLES) {
			ResultSet tables = connection.getMetaData().getTables(null, null, table.toUpperCase(), new String[] { "TABLE" });
			if (!tables.next()) {
				return false;
			}
		}
		return true;
	}

	private void initialise() throws SQLException {

		int tableCount = 0;
		for (String createTable : SQL_CREATE_TABLES) {
			PreparedStatement ps = connection.prepareStatement(createTable);
			ps.execute();
			tableCount++;
		}
		System.out.println("Initialised database schema with " + tableCount + " tables");
	}

	public List<String> listSessionsInDatabase(String username) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_SELECT_SESSIONS_BY_USERNAME);
		ps.setString(1, username);
		ResultSet rs = ps.executeQuery();
		LinkedList<String> sessions = new LinkedList<String>();
		
		while (rs.next()) {
			sessions.add(rs.getString("name"));
		}
		
		return sessions;
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
