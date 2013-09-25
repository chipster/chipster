package fi.csc.microarray.manager;

import it.sauronsoftware.cron4j.Scheduler;

import java.io.File;
import java.io.IOException;
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
import org.eclipse.jetty.security.ConstraintMapping;
import org.eclipse.jetty.security.ConstraintSecurityHandler;
import org.eclipse.jetty.security.HashLoginService;
import org.eclipse.jetty.server.Connector;
import org.eclipse.jetty.server.Handler;
import org.eclipse.jetty.server.handler.DefaultHandler;
import org.eclipse.jetty.server.handler.HandlerCollection;
import org.eclipse.jetty.server.nio.SelectChannelConnector;
import org.eclipse.jetty.util.security.Constraint;
import org.eclipse.jetty.util.security.Password;
import org.eclipse.jetty.util.thread.QueuedThreadPool;
import org.eclipse.jetty.webapp.WebAppContext;
import org.h2.tools.Server;
import org.springframework.core.io.ClassPathResource;
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.jdbc.core.simple.SimpleJdbcInsert;
import org.springframework.jdbc.datasource.DriverManagerDataSource;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingListener;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.MonitoredNodeBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.FeedbackMessage;
import fi.csc.microarray.messaging.message.JobLogMessage;
import fi.csc.microarray.service.KeepAliveShutdownHandler;
import fi.csc.microarray.service.ShutdownCallback;
import fi.csc.microarray.util.Emails;
import fi.csc.microarray.util.MemUtil;

/**
 * Monitoring database and tool for Chipster server system.
 * 
 * Using JdbcTemplate this way for inserts is really slow.
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
    private SimpleJdbcInsert insertAccountTemplate;
    private String feedbackEmail;

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
		");";
	
	/**
	 * Used in admin-web statistics to ignore test accounts 
	 */
	private static final String CREATE_ACCOUNTS_TABLE = 
			"CREATE TABLE IF NOT EXISTS accounts (" +
			"username VARCHAR(200) PRIMARY KEY, " +
			"ignoreInStatistics BOOLEAN DEFAULT FALSE" + 			
			");";
	
	private static final String ADMIN_ROLE = "admin_role";
	
	
	/**
	 * Our route to messaging fabric.
	 */
	private MessagingEndpoint endpoint;
	
	

	/**
	 * 
	 * @throws MicroarrayException
	 * @throws JMSException 
	 * @throws IOException 
	 * @throws Exception 
	 */
	public Manager(String configURL) throws Exception {
		
		// initialise dir and logging
		DirectoryLayout.initialiseServerLayout(Arrays.asList(new String[] {"manager"}), configURL);
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
		logger = Logger.getLogger(Manager.class);
		
		// email for sending feedback
		feedbackEmail = configuration.getString("manager", "admin-email");
		
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
	    this.insertAccountTemplate = new SimpleJdbcInsert(dataSource).withTableName("accounts");

	    // create tables if they do not exist
	    jdbcTemplate.execute(CREATE_JOBS_TABLE);
	    jdbcTemplate.execute(CREATE_ACCOUNTS_TABLE);
		
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
    	
    	// site specific configurations 
    	if (additionalTask != null) {
    		// schedule site specific tasks
    		Scheduler scheduler = new Scheduler();
    		scheduler.schedule("10 0 * * *", additionalTask);
    		scheduler.start();
    	}
    	
    	
		// initialize communications
		this.endpoint = new MessagingEndpoint(this);
		
		// listen for job log messages
		MessagingTopic jobLogTopic = endpoint.createTopic(Topics.Name.JOB_LOG_TOPIC, AccessMode.READ);
		jobLogTopic.setListener(this);
		
	    // listen for feedback messages
        MessagingTopic feedbackTopic = endpoint.createTopic(Topics.Name.AUTHORISED_FEEDBACK_TOPIC, AccessMode.READ);
        feedbackTopic.setListener(this);

		// start h2 web console
		Server h2WebConsoleServer;
		if (startWebConsole) {
			h2WebConsoleServer = Server.createWebServer(new String[] {"-webAllowOthers",  "-webPort", String.valueOf(webConsolePort)});
			h2WebConsoleServer.start();
		}
		
		// start admin
		boolean startAdmin = configuration.getBoolean("manager", "start-admin");
		if (startAdmin) {
    		String s = configuration.getString("manager", "admin-test-account-list").trim();
    		if (!s.isEmpty()) {
    			insertTestAccounts(s.split(","));
    		}
			startAdmin(configuration);
		}
		
		// create keep-alive thread and register shutdown hook
		KeepAliveShutdownHandler.init(this);
		
		logger.error("manager is up and running [" + ApplicationConstants.VERSION + "]");
		logger.info("[mem: " + MemUtil.getMemInfo() + "]");
	}

	private void insertTestAccounts(String[] accounts) {
			
		for (String account : accounts) {
			try {
				
				// add test accounts if they don't exist
				if (jdbcTemplate.queryForList("SELECT * from accounts where username like ?", new String[] {account}).size() == 0) {
					logger.info("adding " + account + " to test accounts");
					
					Map<String, Object> parameters = new HashMap<String, Object>(); 
					parameters.put("username", account);
					parameters.put("ignoreInStatistics", true);

					this.insertAccountTemplate.execute(parameters);
				}
			} catch (Exception e) {
				logger.warn("could not insert test account: " + account, e);
			}
		}
	}

	private void startAdmin(Configuration configuration) throws IOException,
			Exception {
		org.eclipse.jetty.server.Server adminServer = new org.eclipse.jetty.server.Server();
		adminServer.setThreadPool(new QueuedThreadPool());
		Connector connector = new SelectChannelConnector();
		connector.setServer(adminServer);
		connector.setPort(configuration.getInt("manager", "admin-port"));
		adminServer.setConnectors(new Connector[]{ connector });
		
		Constraint constraint = new Constraint();
		constraint.setName(Constraint.__BASIC_AUTH);
		constraint.setRoles(new String[] {ADMIN_ROLE});
		constraint.setAuthenticate(true);
		
		ConstraintMapping cm = new ConstraintMapping();
		cm.setConstraint(constraint);
		cm.setPathSpec("/*");
		
		HashLoginService loginService = new HashLoginService("Please enter Chipster Admin username and password");
		loginService.update(configuration.getString("manager", "admin-username"), 
				new Password(configuration.getString("manager", "admin-password")), 
				new String[] {ADMIN_ROLE});
		
		ConstraintSecurityHandler sh = new ConstraintSecurityHandler();
		sh.setLoginService(loginService);
		sh.addConstraintMapping(cm);
		
		WebAppContext context = new WebAppContext();
		context.setWar(new File(DirectoryLayout.getInstance().getWebappsDir(), "admin-web.war").getAbsolutePath());
        context.setContextPath("/");
		
//        context.setDescriptor(new ClassPathResource("WebContent/WEB-INF/web.xml").getURI().toString());
//        context.setResourceBase(new ClassPathResource("WebContent").getURI().toString());
//        context.setContextPath("/");
//        context.setParentLoaderPriority(true);
				
        context.setHandler(sh);
		HandlerCollection handlers = new HandlerCollection();
		handlers.setHandlers(new Handler[] {context, new DefaultHandler()});
				
		adminServer.setHandler(handlers);
        adminServer.start();
	}

	public String getName() {
		return "manager";
	}


	/**
	 * Process incoming message.  
	 */
	public void onChipsterMessage(ChipsterMessage chipsterMessage) {
		
		if (chipsterMessage instanceof JobLogMessage) {
		    // log information about some job ran by a user
	        JobLogMessage jobLogMessage = (JobLogMessage)chipsterMessage;
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
		} else if (chipsterMessage instanceof FeedbackMessage) {
		    // user gives feedback after seeing an error message
		    FeedbackMessage feedback = (FeedbackMessage) chipsterMessage;
		    logger.info("Feedback received: " + feedback.getDetails());
		    
		    // formulate an email
		    String replyEmail = !feedback.getEmail().equals("") ?
		            feedback.getEmail() : "[not available]";
		    String sessURL = !feedback.getSessionURL().equals("") ? 
		            feedback.getSessionURL() : "[not available]";
		    String emailBody =
		        feedback.getDetails() + "\n\n" +
		        "username: " + feedback.getUsername() + "\n" +
		        "email: " + replyEmail + "\n" +
		        "session file: " + sessURL + "\n\n" +
		        "Download the session file as .zip and open it in Chipster using magic shortcut SHIFT-CTRL-ALT-O\n\n";		    
		    for (String[] log : feedback.getLogs()) {
                emailBody += log[0] + ": " + log[1] + "\n";
            }
		    
		    String subject = "Help request from " + feedback.getUsername();
		    // send the email
		    Emails.sendEmail(feedbackEmail, !feedback.getEmail().equals("") ? feedback.getEmail() : null, subject, emailBody);
		    
		} else {
	        logger.warn("Got other than JobLogMessage: " + chipsterMessage.toString());
	        return; 
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
