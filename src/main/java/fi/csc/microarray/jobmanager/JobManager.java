package fi.csc.microarray.jobmanager;

import java.util.Arrays;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.messaging.JMSMessagingEndpoint;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingListener;
import fi.csc.microarray.messaging.MonitoredNodeBase;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.service.KeepAliveShutdownHandler;
import fi.csc.microarray.service.ShutdownCallback;
import fi.csc.microarray.util.SystemMonitorUtil;

public class JobManager extends MonitoredNodeBase implements MessagingListener, ShutdownCallback {

	private static Logger logger;

	private MessagingEndpoint endpoint;
	

	public JobManager(String configURL) throws Exception {
		
		// initialise dir, config and logging
		DirectoryLayout.initialiseServerLayout(
		        Arrays.asList(new String[] {"jobmanager"}), configURL);
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();

		logger = Logger.getLogger(JobManager.class);
		logger.info("starting jobmanager service...");

		// initialize communications
		this.endpoint = new JMSMessagingEndpoint(this);
		
//		MessagingTopic analyseTopic = endpoint.createTopic(Topics.Name.AUTHORIZED_MANAGED_REQUEST_TOPIC, AccessMode.READ);
//		analyseTopic.setListener(this);

		

		// create keep-alive thread and register shutdown hook
		KeepAliveShutdownHandler.init(this);
		
		logger.info("jobmanager is up and running [" + ApplicationConstants.VERSION + "]");
		logger.info("[mem: " + SystemMonitorUtil.getMemInfo() + "]");
	}
	
	
	@Override
	public void onChipsterMessage(ChipsterMessage msg) {
		// TODO Auto-generated method stub
	}
	
	@Override
	public String getName() {
		return "jobmanager";
	}

	@Override
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
