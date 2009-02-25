package fi.csc.microarray.webstart;

import org.apache.log4j.Logger;
import org.mortbay.jetty.Connector;
import org.mortbay.jetty.Server;
import org.mortbay.jetty.nio.SelectChannelConnector;
import org.mortbay.jetty.servlet.Context;
import org.mortbay.jetty.servlet.DefaultServlet;
import org.mortbay.jetty.servlet.ServletHolder;
import org.mortbay.thread.QueuedThreadPool;

import fi.csc.microarray.ApplicationConstants;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.util.MemUtil;

public class WebstartJettyServer {
	
	private static final int PORT_NUMBER = 8081;
	/**
	 * Logger for this class
	 */
	private static Logger logger;

	static {
		try {
			DirectoryLayout.initialiseServerLayout().getConfiguration();			
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		logger = Logger.getLogger(WebstartJettyServer.class);
	}
	
	private Server jettyInstance;
	
	public static void main(String[] args) throws Exception {

		System.setProperty("DEBUG", "true");
		WebstartJettyServer server = new WebstartJettyServer();
		server.start();
		server.jettyInstance.join();
	}
	
	public void start() {
	
		try {
			if (logger.isDebugEnabled()) {
				System.setProperty("DEBUG", "true");
			}
			jettyInstance = new Server();
			jettyInstance.setThreadPool(new QueuedThreadPool());
			Connector connector = new SelectChannelConnector();
			connector.setServer(jettyInstance);
			connector.setPort(PORT_NUMBER);
			jettyInstance.setConnectors(new Connector[]{ connector });

			Context wsRoot = new Context(jettyInstance, "/", false, false);
			wsRoot.setResourceBase("web-content/");
			wsRoot.addServlet(new ServletHolder(new DefaultServlet()), "/*");

			jettyInstance.start();

			logger.info("webstart is up and running [" + ApplicationConstants.NAMI_VERSION + "]");
			logger.info("[mem: " + MemUtil.getMemInfo() + "]");

		} catch (Exception e) {
			logger.error(e, e);
		}
	}
	
	public boolean isRunning() {
		return jettyInstance.isRunning();
	}
}
