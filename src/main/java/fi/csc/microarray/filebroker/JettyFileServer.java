package fi.csc.microarray.filebroker;

import org.eclipse.jetty.server.Connector;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.nio.SelectChannelConnector;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;
import org.eclipse.jetty.util.thread.QueuedThreadPool;

import fi.csc.microarray.config.DirectoryLayout;

public class JettyFileServer {

	private Server jettyInstance;
	private AuthorisedUrlRepository urlRepository;
	
	public JettyFileServer(AuthorisedUrlRepository urlRepository) {
		this.urlRepository = urlRepository;
	}
	
	public void start(String resourceBase, int port) throws Exception {
		
		if (DirectoryLayout.getInstance().getConfiguration().getBoolean("filebroker", "jetty-debug")) {
			System.setProperty("DEBUG", "true");
		}
		
		jettyInstance = new Server();
		jettyInstance.setThreadPool(new QueuedThreadPool());
		Connector connector = new SelectChannelConnector();
		connector.setServer(jettyInstance);
		connector.setPort(port);
		jettyInstance.setConnectors(new Connector[]{ connector });

		ServletContextHandler root = new ServletContextHandler(jettyInstance, "/", false, false);
		root.getInitParams().put("org.eclipse.jetty.servlet.Default.aliases", "true");
		root.setResourceBase(resourceBase);
		root.addServlet(new ServletHolder(new RestServlet(urlRepository, urlRepository.getRootUrl())), "/*");
		jettyInstance.start();
	}
	
	public boolean isRunning() {
		return jettyInstance.isRunning();
	}
}
