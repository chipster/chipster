package fi.csc.microarray.webstart;

import java.io.File;
import java.util.Arrays;

import org.apache.log4j.Logger;
import org.eclipse.jetty.server.Connector;
import org.eclipse.jetty.server.Handler;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.handler.HandlerCollection;
import org.eclipse.jetty.server.nio.SelectChannelConnector;
import org.eclipse.jetty.servlet.DefaultServlet;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;
import org.eclipse.jetty.util.thread.QueuedThreadPool;
import org.eclipse.jetty.webapp.WebAppContext;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.service.KeepAliveShutdownHandler;
import fi.csc.microarray.service.ShutdownCallback;
import fi.csc.microarray.util.MemUtil;

public class WebstartJettyServer implements ShutdownCallback {
	
	/**
	 * Logger for this class
	 */
	private static Logger logger;
	
	private Server jettyInstance;
	
	public static void main(String[] args) throws Exception {

		System.setProperty("DEBUG", "true");
		WebstartJettyServer server = new WebstartJettyServer();
		server.start();
		server.jettyInstance.join();
	}
	
	public void start() {
	
		try {
			// initialise config and logging
			DirectoryLayout.initialiseServerLayout(Arrays.asList(new String[] {"webstart"}));
			Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
			logger = Logger.getLogger(WebstartJettyServer.class);
			if (logger.isDebugEnabled()) {
				System.setProperty("DEBUG", "true");
			}

			// initialise jetty
			jettyInstance = new Server();
			jettyInstance.setThreadPool(new QueuedThreadPool());
			Connector connector = new SelectChannelConnector();
			connector.setServer(jettyInstance);
			connector.setPort(configuration.getInt("webstart", "port"));
			jettyInstance.setConnectors(new Connector[]{ connector });

			// webstart front page
			ServletContextHandler wsRoot = new ServletContextHandler(jettyInstance, "/", false, false);
			wsRoot.setResourceBase(DirectoryLayout.WEB_ROOT + "/");
			wsRoot.addServlet(new ServletHolder(new DefaultServlet()), "/*");

			// tooleditor web app
			WebAppContext context = new WebAppContext();
			context.setWar(new File(DirectoryLayout.getInstance().getWebappsDir(), "tool-editor.war").getAbsolutePath());
	        context.setContextPath("/tool-editor");
			
	        // seems to work like this, no idea if this is correct way to do it
			HandlerCollection handlers = new HandlerCollection();
			handlers.setHandlers(new Handler[] {context, wsRoot});
			jettyInstance.setHandler(handlers);
			
			jettyInstance.start();

			// create keep-alive thread and register shutdown hook
			KeepAliveShutdownHandler.init(this);
			
			logger.info("webstart is up and running [" + ApplicationConstants.VERSION + "]");
			logger.info("[mem: " + MemUtil.getMemInfo() + "]");

		} catch (Exception e) {
			e.printStackTrace();
			logger.error(e, e);
		}
	}
	
	public boolean isRunning() {
		return jettyInstance.isRunning();
	}

	public void shutdown() {
		logger.info("shutdown requested");
		logger.info("shutting down");
	}


}
