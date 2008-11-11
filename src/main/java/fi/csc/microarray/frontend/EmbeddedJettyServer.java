package fi.csc.microarray.frontend;

import org.mortbay.jetty.Connector;
import org.mortbay.jetty.Server;
import org.mortbay.jetty.nio.SelectChannelConnector;
import org.mortbay.jetty.servlet.Context;
import org.mortbay.jetty.servlet.ServletHolder;
import org.mortbay.thread.QueuedThreadPool;

import fi.csc.microarray.MicroarrayConfiguration;
import fi.csc.microarray.util.rest.RestServlet;

public class EmbeddedJettyServer {

	static {
		try {
			MicroarrayConfiguration.loadConfiguration();
		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}
	
	private Server jettyInstance;
	
	public static void main(String[] args) throws Exception {

		System.setProperty("DEBUG", "true");
		EmbeddedJettyServer server = new EmbeddedJettyServer();
		server.start("/home/akallio/eclipse-workspace/microarray/nami-work-files/file-repository/", "/fileserver", 8080);
		server.jettyInstance.join();
	}
	
	public void start(String resourceBase, String contextPath, int port) throws Exception {
		
		if ("true".equals(MicroarrayConfiguration.getValue("filebroker", "jettyDebug"))) {
			System.setProperty("DEBUG", "true");
		}
		
		jettyInstance = new Server();
		jettyInstance.setThreadPool(new QueuedThreadPool());
		Connector connector = new SelectChannelConnector();
		connector.setServer(jettyInstance);
		connector.setPort(port);
		jettyInstance.setConnectors(new Connector[]{ connector });

		Context root = new Context(jettyInstance, contextPath, false, false);
		root.setResourceBase(resourceBase);
		root.addServlet(new ServletHolder(new RestServlet()), "/*");
		jettyInstance.start();
	}
	
	public boolean isRunning() {
		return jettyInstance.isRunning();
	}
}
