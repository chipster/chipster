package fi.csc.microarray.filebroker;

import org.eclipse.jetty.server.Connector;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.nio.SelectChannelConnector;
import org.eclipse.jetty.server.ssl.SslSelectChannelConnector;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;
import org.eclipse.jetty.util.thread.QueuedThreadPool;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.util.KeyAndTrustManager;

public class JettyFileServer {

	private Server jettyInstance;
	private AuthorisedUrlRepository urlRepository;
	private DerbyMetadataServer metadataServer;
	
	public JettyFileServer(AuthorisedUrlRepository urlRepository, DerbyMetadataServer metadataServer) {
		this.urlRepository = urlRepository;
		this.metadataServer = metadataServer;
	}
	
	public void start(String resourceBase, int port, String protocol) throws Exception {
		
		if (DirectoryLayout.getInstance().getConfiguration().getBoolean("filebroker", "jetty-debug")) {
			System.setProperty("org.eclipse.jetty.LEVEL", "DEBUG");
		}
		
		jettyInstance = new Server();
		jettyInstance.setThreadPool(new QueuedThreadPool());
		
		Connector connector;
		switch (protocol) {
		case "http":
			connector= new SelectChannelConnector();
			break;
			
		case "https":
			Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
			connector = new SslSelectChannelConnector(KeyAndTrustManager.createSslContextFactory(
					configuration.getString("security", "keystore"),
					configuration.getString("security", "keypass"), 
					configuration.getString("security", "keyalias"), 
					configuration.getString("security", "master-keystore")
			));
			break;
			
		default:
			throw new IllegalArgumentException("unsupported protocol: " + protocol + " (supported are http and https)");
		}
		connector.setServer(jettyInstance);
		connector.setPort(port);
		jettyInstance.setConnectors(new Connector[]{ connector });

		ServletContextHandler root = new ServletContextHandler(jettyInstance, "/", false, false);
		root.getInitParams().put("org.eclipse.jetty.servlet.Default.aliases", "true");
		root.setResourceBase(resourceBase);
		root.addServlet(new ServletHolder(new RestServlet(urlRepository.getRootUrl(), urlRepository, metadataServer)), "/*");
		jettyInstance.start();
	}
	
	public boolean isRunning() {
		return jettyInstance.isRunning();
	}
}
