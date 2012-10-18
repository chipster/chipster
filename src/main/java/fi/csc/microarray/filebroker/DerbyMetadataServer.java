package fi.csc.microarray.filebroker;

import java.net.URL;
import java.sql.Connection;
import java.sql.DriverManager;
import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.util.IOUtils;

public class DerbyMetadataServer {

	private Connection connection = null;

	private static String sessionDbTable = "sessions";
	
	private static String SQL_CREATE_TABLE = "CREATE TABLE " + sessionDbTable + " (" +
			"id INTEGER NOT NULL GENERATED ALWAYS AS IDENTITY (START WITH 1, INCREMENT BY 1) PRIMARY KEY,  " +
			"name VARCHAR(200),  " +
			"username VARCHAR(200), " +
			"url VARCHAR(200))";

	private static String SQL_SELECT_SESSIONS_BY_USERNAME = "SELECT * FROM " + sessionDbTable + " WHERE username = ?";
	private static String SQL_INSERT_SESSION  = "INSERT INTO " + sessionDbTable + " (name, username, url) VALUES (?, ?, ?)";
	private static String SQL_DELETE_SESSION  = "DELETE FROM " + sessionDbTable + " WHERE name = ? AND username = ?";
	
	
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
		ResultSet tables = connection.getMetaData().getTables(null, null, sessionDbTable.toUpperCase(), new String[] { "TABLE" });
		return tables.next();
	}

	private void initialise() throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_CREATE_TABLE);
        ps.execute();
        System.out.println("Table " + sessionDbTable + " initialised");
	}

	public List<String> listSessionsInDatabase(String username) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_SELECT_SESSIONS_BY_USERNAME);
		ResultSet rs = ps.executeQuery();
		LinkedList<String> sessions = new LinkedList<String>();
		
		while (rs.next()) {
			sessions.add(rs.getString("name"));
		}
		
		return sessions;
	}

	public void addSessionToDatabase(String username, String name, URL url) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_INSERT_SESSION);
		ps.setString(1, name);
		ps.setString(2, username);
		ps.setString(3, IOUtils.getFilenameWithoutPath(url));
		ps.execute();
	}

	public void removeSessionFromDatabase(String username, String name) throws SQLException {
		PreparedStatement ps = connection.prepareStatement(SQL_DELETE_SESSION);
		ps.setString(1, name);
		ps.setString(2, username);
		ps.execute();
	}
	
	
}
