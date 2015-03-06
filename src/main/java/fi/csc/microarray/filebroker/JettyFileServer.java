package fi.csc.microarray.filebroker;

import org.eclipse.jetty.server.Connector;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.nio.SelectChannelConnector;
import org.eclipse.jetty.server.ssl.SslSocketConnector;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;
import org.eclipse.jetty.util.ssl.SslContextFactory;
import org.eclipse.jetty.util.thread.QueuedThreadPool;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.util.KeyAndTrustManager;

public class JettyFileServer {

	private Server jettyInstance;
	private AuthorisedUrlRepository urlRepository;
	private DerbyMetadataServer metadataServer;
	private DiskCleanUp cacheCleanUp;
	
	public JettyFileServer(AuthorisedUrlRepository urlRepository, DerbyMetadataServer metadataServer, DiskCleanUp cacheCleanUp) {
		this.urlRepository = urlRepository;
		this.metadataServer = metadataServer;
		this.cacheCleanUp = cacheCleanUp;
	}

	public void start(String resourceBase, int port, String protocol) throws Exception {
		
		if (DirectoryLayout.getInstance().getConfiguration().getBoolean("filebroker", "jetty-debug")) {
			System.setProperty("org.eclipse.jetty.LEVEL", "DEBUG");
		}
		
		jettyInstance = new Server();
		jettyInstance.setThreadPool(new QueuedThreadPool());
		
		Connector connector;
		SslContextFactory contextFactory = null;
		
		switch (protocol) {
		case "http":
			connector= new SelectChannelConnector();
			break;
			
		case "https":
			Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
			
			String[] protocols = configuration.getString("filebroker", "ssl-protocol-version").split(",");
			
			contextFactory = KeyAndTrustManager.createSslContextFactory(
					configuration.getString("security", "filebroker-keystore"),
					configuration.getString("security", "storepass"),
					protocols
			);
			connector = new SslSocketConnector(contextFactory);			
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
		root.addServlet(new ServletHolder(new RestServlet(urlRepository.getRootUrl(), urlRepository, metadataServer, cacheCleanUp)), "/*");
		jettyInstance.start();
	}
	
	public boolean isRunning() {
		return jettyInstance.isRunning();
	}
}
