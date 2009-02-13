package fi.csc.microarray.filebroker;

import java.util.TimerTask;

import org.apache.log4j.Logger;

public class JettyCheckTimerTask extends TimerTask {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(JettyCheckTimerTask.class);

	private EmbeddedJettyServer jettyServer;

	public JettyCheckTimerTask(EmbeddedJettyServer jettyInstance) {
		this.jettyServer = jettyInstance; 		
	}
	
	@Override
	public void run() {
		if (!jettyServer.isRunning()) {
			logger.error("jetty server is not running");
		}

	}

}
