package fi.csc.microarray.manager;

import java.io.IOException;
import java.sql.SQLException;
import java.util.HashMap;
import java.util.Map;

import javax.jms.JMSException;

import org.apache.log4j.Logger;
import org.h2.tools.Server;
import org.springframework.jdbc.core.JdbcTemplate;
import org.springframework.jdbc.core.simple.SimpleJdbcInsert;
import org.springframework.jdbc.datasource.DriverManagerDataSource;

import fi.csc.microarray.ApplicationConstants;
import fi.csc.microarray.MicroarrayConfiguration;
import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingListener;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.MonitoredNodeBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.message.JobLogMessage;
import fi.csc.microarray.messaging.message.NamiMessage;
import fi.csc.microarray.util.MemUtil;
import fi.csc.microarray.util.config.ConfigurationLoader.OldConfigurationFormatException;

/**
 * 
 * @author hupponen
 */
public class Manager extends MonitoredNodeBase implements MessagingListener {
	

	
	
	/**
	 * Loggers.
	 */
	private final Logger logger;

    private JdbcTemplate jdbcTemplate;
    private SimpleJdbcInsert insertJobTemplate;

	private static final String DATABASE = "database/chipster-manager";
	private static final boolean START_WEB_CONSOLE = true;
	private static final String DB_USERNAME = "chipster";
	private static final String DB_PASSWORD = "";
    
    
	// TODO index, unique keys
	private static final String CREATE_JOBS_TABLE = 
		"CREATE TABLE IF NOT EXISTS jobs (" +
			"jobId VARCHAR(100), " + 
			"operation VARCHAR(100), " +
			"status VARCHAR(20), " + 
			"starttime DATETIME DEFAULT NULL, " + 
			"endtime DATETIME DEFAULT NULL, " +
			"wallclockTime INT DEFAULT NULL, " +
			"errorMessage TEXT DEFAULT NULL, " +
			"outputText TEXT DEFAULT NULL, " + 
			"username VARCHAR(20), " +
			"compHost VARCHAR(50)" +
			");";
	
	
	
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
	 * @throws OldConfigurationFormatException 
	 * @throws IOException 
	 * @throws ClassNotFoundException 
	 * @throws SQLException 
	 */
	public Manager() throws MicroarrayException, JMSException, IOException, OldConfigurationFormatException, ClassNotFoundException, SQLException {
		
		MicroarrayConfiguration.loadConfiguration();
		logger = Logger.getLogger(Manager.class);
		logger.info("Starting manager...");
		
		
		// initialize database connection
		DriverManagerDataSource dataSource = new DriverManagerDataSource();
		dataSource.setDriverClassName("org.h2.Driver");
		dataSource.setUrl("jdbc:h2:" + DATABASE);
		dataSource.setUsername(DB_USERNAME);
		dataSource.setPassword(DB_PASSWORD);
		
        this.jdbcTemplate = new JdbcTemplate(dataSource);
	    this.insertJobTemplate = new SimpleJdbcInsert(dataSource).withTableName("jobs");

	    // create tables if they do not exist
	    jdbcTemplate.execute(CREATE_JOBS_TABLE);
		
		// initialize communications
		this.endpoint = new MessagingEndpoint(this);
		
		MessagingTopic managerTopic = endpoint.createTopic(Topics.Name.MANAGER_TOPIC, AccessMode.READ);
		managerTopic.setListener(this);

		// start web console
		Server server;
		if (START_WEB_CONSOLE) {
			server = Server.createWebServer(new String[] {"-webAllowOthers"});
			server.start();
		}
		
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
		
	    Map<String, Object> parameters = new HashMap<String, Object>();
	    parameters.put("jobId", jobLogMessage.getJobId());
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
	
	}


}
