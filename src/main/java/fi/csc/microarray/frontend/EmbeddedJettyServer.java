package fi.csc.microarray.frontend;

import org.apache.log4j.Logger;
import org.mortbay.jetty.Connector;
import org.mortbay.jetty.Server;
import org.mortbay.jetty.bio.SocketConnector;
import org.mortbay.jetty.servlet.Context;
import org.mortbay.jetty.servlet.ServletHolder;

import fi.csc.microarray.MicroarrayConfiguration;
import fi.csc.microarray.util.rest.RestServlet;

public class EmbeddedJettyServer {
	/**
	 * Logger for this class
	 */
	private static Logger logger;

	static {
		try {
			MicroarrayConfiguration.loadConfiguration();
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
		logger = Logger.getLogger(EmbeddedJettyServer.class);
	}
	
	private Server jettyInstance;
	
	public static void main(String[] args) throws Exception {

		System.setProperty("DEBUG", "true");
		EmbeddedJettyServer server = new EmbeddedJettyServer();
		server.start("/home/akallio/eclipse-workspace/microarray/nami-work-files/file-repository/", "/fileserver", 8080);
		server.jettyInstance.join();
	}
	
	public void start(String resourceBase, String contextPath, int port) throws Exception {
		
		if (logger.isDebugEnabled()) {
			System.setProperty("DEBUG", "true");
		}
		jettyInstance = new Server();
		Connector connector = new SocketConnector();
		connector.setServer(jettyInstance);
		connector.setPort(port);
		jettyInstance.setConnectors(new Connector[]{connector});

		Context root = new Context(jettyInstance, contextPath, false, false);
		root.setResourceBase(resourceBase);
		root.addServlet(new ServletHolder(new RestServlet()), "/*");
		// TODO verify if FilenameGuardFilter is causing trouble ("dispatch failed" in some cases)
		//root.addFilter(new FilterHolder(new FilenameGuardFilter()), "/*", 1); 
		jettyInstance.start();
	}
	
	public boolean isRunning() {
		return jettyInstance.isRunning();
	}
}
