package fi.csc.microarray.filebroker;

import org.mortbay.jetty.Connector;
import org.mortbay.jetty.Server;
import org.mortbay.jetty.nio.SelectChannelConnector;
import org.mortbay.jetty.servlet.Context;
import org.mortbay.jetty.servlet.ServletHolder;
import org.mortbay.thread.QueuedThreadPool;

import fi.csc.microarray.config.DirectoryLayout;

public class JettyFileServer {

	private Server jettyInstance;
	private AuthorisedUrlRepository urlRepository;
	
	public JettyFileServer(AuthorisedUrlRepository urlRepository) {
		this.urlRepository = urlRepository;
	}
	
	@SuppressWarnings("unchecked")
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

		Context root = new Context(jettyInstance, "/", false, false);
		root.getInitParams().put("org.mortbay.jetty.servlet.Default.aliases", "true");
		root.setResourceBase(resourceBase);
		root.addServlet(new ServletHolder(new RestServlet(urlRepository, urlRepository.getRootUrl())), "/*");
		jettyInstance.start();
	}
	
	public boolean isRunning() {
		return jettyInstance.isRunning();
	}
}
