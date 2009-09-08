package fi.csc.microarray.manager;

import it.sauronsoftware.cron4j.Scheduler;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Arrays;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.Map;
import java.util.Timer;
import java.util.TimerTask;

import javax.jms.JMSException;

import org.apache.log4j.Logger;
import org.h2.tools.Server;
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.jdbc.core.simple.SimpleJdbcInsert;
import org.springframework.jdbc.datasource.DriverManagerDataSource;

import fi.csc.microarray.ApplicationConstants;
import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingListener;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.MonitoredNodeBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.message.JobLogMessage;
import fi.csc.microarray.messaging.message.NamiMessage;
import fi.csc.microarray.service.KeepAliveShutdownHandler;
import fi.csc.microarray.service.ShutdownCallback;
import fi.csc.microarray.util.MemUtil;

/**
 * Monitoring database and tool for Chipster server system.
 * 
 * @author Taavi Hupponen
 */
public class Manager extends MonitoredNodeBase implements MessagingListener, ShutdownCallback {
	

	private class BackupTimerTask extends TimerTask {

		private File baseBackupDir;
		
		public BackupTimerTask(File backupDir) {
			this.baseBackupDir = backupDir;
		}
		
		
		@Override
		public void run() {
			logger.info("Creating database backup");
			DateFormat df = new SimpleDateFormat("yyyy-MM-dd_mm:ss.SSS");
			String fileName = baseBackupDir.getAbsolutePath() + File.separator + "chipster-manager-db-backup-" + df.format(new Date());
			String sql = "SCRIPT TO '" + fileName + ".zip' COMPRESSION ZIP";
			jdbcTemplate.execute(sql);
		}
		
	}
	
	
	/**
	 * Loggers.
	 */
	private final Logger logger;

    private JdbcTemplate jdbcTemplate;
    private SimpleJdbcInsert insertJobTemplate;

	// TODO index, unique keys
	private static final String CREATE_JOBS_TABLE = 
		"CREATE TABLE IF NOT EXISTS jobs (" +
		"id VARCHAR(100) PRIMARY KEY, " + 
		"operation VARCHAR(100), " +
		"status VARCHAR(100), " + 
		"starttime DATETIME DEFAULT NULL, " + 
		"endtime DATETIME DEFAULT NULL, " +
		"wallclockTime INT DEFAULT NULL, " +
		"errorMessage TEXT DEFAULT NULL, " +
		"outputText TEXT DEFAULT NULL, " + 
		"username VARCHAR(200), " +
		"compHost VARCHAR(200)" +
		"); ";
		//"CREATE UNIQUE INDEX IF NOT EXISTS jobIdIndex on jobs(id); ";
	
	
	
	/**
	 * Our route to messaging fabric.
	 */
	private MessagingEndpoint endpoint;
	
	

	/**
	 * 
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
	public Manager() throws MicroarrayException, JMSException, IOException, IllegalConfigurationException, ClassNotFoundException, SQLException {
		
		// initialise dir and logging
		DirectoryLayout.initialiseServerLayout(Arrays.asList(new String[] {"manager"}));
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
		logger = Logger.getLogger(Manager.class);
		
		// initialize database connection
		logger.info("starting manager...");
		String dbDriver = configuration.getString("manager", "jdbc-driver");
		String dbUrl = configuration.getString("manager", "database-url");
		boolean startWebConsole = configuration.getBoolean("manager", "start-web-console");
		String dbUsername = configuration.getString("manager", "database-username");
	    String dbPassword = configuration.getString("manager", "database-password");
	    int webConsolePort = configuration.getInt("manager", "web-console-port");

		
		// FIXME we should retrieve database directory from DirectoryLayout
		DriverManagerDataSource dataSource = new DriverManagerDataSource();
		dataSource.setDriverClassName(dbDriver);
		dataSource.setUrl(dbUrl);
		dataSource.setUsername(dbUsername);
		dataSource.setPassword(dbPassword);
		
        this.jdbcTemplate = new JdbcTemplate(dataSource);
	    this.insertJobTemplate = new SimpleJdbcInsert(dataSource).withTableName("jobs");

	    // create tables if they do not exist
	    jdbcTemplate.execute(CREATE_JOBS_TABLE);
		
	    // schedule backups
	    int backupInterval = configuration.getInt("manager", "backup-interval");
	    String backupTimeString = configuration.getString("manager", "backup-time");
	    int startHour = Integer.parseInt(backupTimeString.split(":")[0]);
	    int startMinute = Integer.parseInt(backupTimeString.split(":")[1]);
	    Calendar firstBackupTime = Calendar.getInstance();
	    if (firstBackupTime.get(Calendar.HOUR_OF_DAY) > startHour || 
	    		(firstBackupTime.get(Calendar.HOUR_OF_DAY) == startHour && 
	    				firstBackupTime.get(Calendar.MINUTE) >= startMinute)) {
	    	firstBackupTime.add(Calendar.DATE, 1);
	    }
    	firstBackupTime.set(Calendar.HOUR_OF_DAY, startHour);
    	firstBackupTime.set(Calendar.MINUTE, startMinute);
    	firstBackupTime.set(Calendar.SECOND, 0);
    	firstBackupTime.set(Calendar.MILLISECOND, 0);
    	logger.info("Next database backup is scheduled at " + firstBackupTime.getTime().toString());
    	
	    Timer timer = new Timer("chipster-manager-backup", true);
	    File backupDir = DirectoryLayout.getInstance().getBackupDir();
    	timer.scheduleAtFixedRate(new BackupTimerTask(backupDir), firstBackupTime.getTime(), backupInterval*60*60*1000);
	    
    	
	    // schedule additional tasks
    	// TODO add to configs
    	TimerTask additionalTask = null;
    	try {
    		additionalTask = (TimerTask) Class.forName("fi.csc.chipster.manager.AskareLogTimerTask").getConstructor(JdbcTemplate.class).newInstance(jdbcTemplate);
    	} catch (Exception e) {
    		logger.info("could not load additional tasks");
    	}
    	if (additionalTask != null) {
    		Scheduler scheduler = new Scheduler();
    		scheduler.schedule("10 0 * * *", additionalTask);
    		scheduler.start();
    	}
    	
    	
		// initialize communications
		this.endpoint = new MessagingEndpoint(this);
		
		MessagingTopic managerTopic = endpoint.createTopic(Topics.Name.MANAGER_TOPIC, AccessMode.READ);
		managerTopic.setListener(this);

		// start web console
		Server server;
		if (startWebConsole) {
			server = Server.createWebServer(new String[] {"-webAllowOthers",  "-webPort", String.valueOf(webConsolePort)});
			server.start();
		}
		
		// create keep-alive thread and register shutdown hook
		KeepAliveShutdownHandler.init(this);
		
		logger.error("manager is up and running [" + ApplicationConstants.NAMI_VERSION + "]");
		logger.info("[mem: " + MemUtil.getMemInfo() + "]");
	}
	

	public String getName() {
		return "manager";
	}


	/**
	 * Process incoming message.  
	 */
	public void onNamiMessage(NamiMessage namiMessage) {
		
		if (!(namiMessage instanceof JobLogMessage)) {
			logger.warn("Got other than JobLogMessage: " + namiMessage.toString());
			return;
		}
		
		JobLogMessage jobLogMessage = (JobLogMessage)namiMessage;
		try {
		    Map<String, Object> parameters = new HashMap<String, Object>();
		    parameters.put("id", jobLogMessage.getJobId());
		    parameters.put("operation", jobLogMessage.getOperation());
			parameters.put("status", jobLogMessage.getState().toString()); 
			parameters.put("starttime", jobLogMessage.getStartTime()); 
			parameters.put("endtime", jobLogMessage.getEndTime());
			parameters.put("wallclockTime", (jobLogMessage.getEndTime().getTime() - jobLogMessage.getStartTime().getTime()) / 1000);
			parameters.put("errorMessage", jobLogMessage.getErrorMessage());
			parameters.put("outputText", jobLogMessage.getOutputText()); 
			parameters.put("username", jobLogMessage.getUsername());
			parameters.put("compHost", jobLogMessage.getCompHost());
			
			this.insertJobTemplate.execute(parameters);
		} catch (Exception e) {
			logger.error("Could not insert log entry", e);
		}
	}

	public void shutdown() {
		logger.info("shutdown requested");

		// close messaging endpoint
		try {
			this.endpoint.close();
		} catch (JMSException e) {
			logger.error("closing messaging endpoint failed", e);
		}

		logger.info("shutting down");
	}
	
}
