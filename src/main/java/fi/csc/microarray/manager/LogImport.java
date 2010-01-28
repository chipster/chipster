package fi.csc.microarray.manager;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.text.DateFormat;
import java.text.ParseException;
import java.text.SimpleDateFormat;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Map;
import java.util.UUID;

import javax.jms.JMSException;

import org.h2.tools.Server;
import org.springframework.dao.DataIntegrityViolationException;
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.jdbc.core.simple.SimpleJdbcInsert;
import org.springframework.jdbc.datasource.DriverManagerDataSource;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;


public class LogImport {

		
		
		/**
		 * Loggers.
		 */

	    private JdbcTemplate jdbcTemplate;
	    private SimpleJdbcInsert insertJobTemplate;

		// TODO index, unique keys
		private static final String CREATE_JOBS_TABLE = 
			"CREATE TABLE IF NOT EXISTS jobs (" +
				"id VARCHAR(100) PRIMARY KEY, " + 
				"operation VARCHAR(200), " +
				"status VARCHAR(200), " + 
				"starttime DATETIME DEFAULT NULL, " + 
				"endtime DATETIME DEFAULT NULL, " +
				"wallclockTime INT DEFAULT NULL, " +
				"errorMessage TEXT DEFAULT NULL, " +
				"outputText TEXT DEFAULT NULL, " + 
				"username VARCHAR(200), " +
				"compHost VARCHAR(500)" +
				"); "; //+
//				"CREATE UNIQUE INDEX IF NOT EXISTS jobIdIndex on jobs(id); ";
		
		
		

		/**
		 * 
		 * @throws SQLException 
		 * @throws MicroarrayException 
		 * @throws JMSException
		 * @throws IOException if creation of working directory fails.
		 * @throws MicroarrayException
		 * @throws JMSException 
		 * @throws IllegalConfigurationException 
		 * @throws IOException 
		 * @throws ClassNotFoundException 
		 * @throws SQLException 
		 */
		public LogImport() throws SQLException {
			
			
			// initialize database connection
//			String dbDriver = MicroarrayConfiguration.getValue("manager", "jdbcDriver");
//			String dbUrl = MicroarrayConfiguration.getValue("manager", "databaseUrl");
//			boolean startWebConsole = "true".equals(MicroarrayConfiguration.getValue("manager", "startWebConsole"));
//			String dbUsername = MicroarrayConfiguration.getValue("manager", "databaseUsername");
//		    String dbPassword = MicroarrayConfiguration.getValue("manager", "databasePassword");
//		    int webConsolePort = Integer.parseInt(MicroarrayConfiguration.getValue("manager", "webConsolePort"));

//			Server tcpServer;
//			tcpServer = Server.createTcpServer();
//			tcpServer.start();
			
			String dbDriver = "org.h2.Driver";
			String dbUrl = "jdbc:h2:~/database/chipster-manager";
			boolean startWebConsole = true;
			String dbUsername = "chipster";
		    String dbPassword = "";
		    int webConsolePort = 8082;

		    
			DriverManagerDataSource dataSource = new DriverManagerDataSource();
			dataSource.setDriverClassName(dbDriver);
			dataSource.setUrl(dbUrl);
			dataSource.setUsername(dbUsername);
			dataSource.setPassword(dbPassword);
			
	        this.jdbcTemplate = new JdbcTemplate(dataSource);
		    this.insertJobTemplate = new SimpleJdbcInsert(dataSource).withTableName("jobs");

		    // create tables if they do not exist
		    jdbcTemplate.execute(CREATE_JOBS_TABLE);
			
			// start web console
			Server server;
			if (startWebConsole) {
				server = Server.createWebServer(new String[] {"-webAllowOthers",  "-webPort", String.valueOf(webConsolePort)});
				server.start();
			}
		}

		public void importLog(File logFile) throws IOException, ParseException {
			DateFormat df = new SimpleDateFormat("EEE MMM dd HH:mm:ss z yyyy");
			LogParser parser = new LogParser();
			Iterable<HashMap<String, String>> entries = parser.parse(logFile);
			int success = 0;
			int integrityFailed = 0;
			int usernameFailed = 0;
			int endtimeFailed = 0;
			int failed = 0;
			for (HashMap<String, String> entry: entries) {
			    try {
			    	// some sanity checks
			    	Date startTime = df.parse(entry.get("starttime"));
			    	Date endTime = df.parse(entry.get("endtime"));
			    	
			    	if (entry.get("username").equals("null")) {
			    		usernameFailed++;
			    		continue;
			    	}
			    	
					Map<String, Object> parameters = new HashMap<String, Object>();
				    parameters.put("id", UUID.randomUUID().toString());
				    parameters.put("operation", entry.get("operation"));
					parameters.put("status", entry.get("status")); 
					parameters.put("starttime", startTime); 
					parameters.put("endtime", endTime);
					parameters.put("wallclockTime", (endTime.getTime() - startTime.getTime()) / 1000);
					parameters.put("errorMessage", null);
					parameters.put("outputText", null); 
					parameters.put("username", entry.get("username"));
					parameters.put("compHost", entry.get("compHost"));
					
					this.insertJobTemplate.execute(parameters);
					success++;
			    } catch (DataIntegrityViolationException dive) {
			    	//System.out.println(dive.toString());
			    	integrityFailed++;
			    }
			    
			    catch (Exception e) {
			    	System.out.println(e.toString());
			    	failed++;
			    	
			    }
			}
			System.out.println("success: " + success + 
					", integrity failed: " + integrityFailed + 
					", username failed: " + usernameFailed +
					", endtime failed: " + endtimeFailed +
					", other failed: " + failed);
		}

		
		

	public void importOldLog(File logFile, String compHost) throws IOException, ParseException {
		DateFormat df = new SimpleDateFormat("E MMM d H:m:s z yyyy");
		LogParser parser = new LogParser();
		Iterable<HashMap<String, String>> entries = parser.parse(logFile);
		int success = 0;
		int integrityFailed = 0;
		int usernameFailed = 0;
		int endtimeFailed = 0;
		int failed = 0;
		for (HashMap<String, String> entry: entries) {
		    try {
		    	// some sanity checks
		    	Date startTime = df.parse(entry.get("starttime"));
		    	Date endTime;
		    	if (!entry.get("endtime").equals("null")) {
		    		endTime = df.parse(entry.get("endtime"));
		    	} else {
		    		// if endtime missing set duration to 5 seconds
		    		Calendar calendar = new GregorianCalendar();
		    		calendar.setTime(startTime);
		    		calendar.add(GregorianCalendar.SECOND, 5);
		    		endTime = calendar.getTime();
		    		endtimeFailed++;
		    	}
		    	
		    	if (entry.get("username").equals("null")) {
		    		usernameFailed++;
		    		continue;
		    	}
		    	
				Map<String, Object> parameters = new HashMap<String, Object>();
			    parameters.put("id", UUID.randomUUID().toString());
			    parameters.put("operation", entry.get("operation"));
				parameters.put("status", entry.get("status")); 
				parameters.put("starttime", startTime); 
				parameters.put("endtime", endTime);
				parameters.put("wallclockTime", (endTime.getTime() - startTime.getTime()) / 1000);
				parameters.put("errorMessage", null);
				parameters.put("outputText", null); 
				parameters.put("username", entry.get("username"));
				parameters.put("compHost", compHost);
				
				this.insertJobTemplate.execute(parameters);
				success++;
		    } catch (DataIntegrityViolationException dive) {
		    	//System.out.println(dive.toString());
		    	integrityFailed++;
		    }
		    
		    catch (Exception e) {
		    	System.out.println(e.toString());
		    	failed++;
		    	
		    }
		}
		System.out.println("success: " + success + 
				", integrity failed: " + integrityFailed + 
				", username failed: " + usernameFailed +
				", endtime failed: " + endtimeFailed +
				", other failed: " + failed);
	}
	
	
	
	public static void main(String[] args) throws SQLException, IOException, IllegalConfigurationException, ParseException {
//		MicroarrayConfiguration.loadConfiguration();
		LogImport logImport = new LogImport();
		logImport.importLog(new File("logfile.log"));
		
	}

}
