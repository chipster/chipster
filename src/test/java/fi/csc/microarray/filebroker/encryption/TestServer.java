package fi.csc.microarray.filebroker.encryption;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.eclipse.jetty.server.Connector;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.ServerConnector;
import org.eclipse.jetty.servlet.DefaultServlet;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;
import org.eclipse.jetty.util.IO;
import org.eclipse.jetty.util.URIUtil;
import org.eclipse.jetty.util.ssl.SslContextFactory;

import sun.net.www.protocol.http.HttpURLConnection;

public class TestServer {

	
	public static void main(String[] args) throws Exception {
		Server jettyInstance = new Server();
		ServerConnector connector = new ServerConnector(jettyInstance);
		connector.setPort(9080);
		SslContextFactory sslContext = new SslContextFactory("keystore.ks");
		sslContext.setKeyStorePassword("microarray");
		ServerConnector sslConnector = new ServerConnector(jettyInstance);
		sslConnector.setPort(9443);
		jettyInstance.setConnectors(new Connector[]{ connector, sslConnector });

		ServletContextHandler root = new ServletContextHandler(jettyInstance, "/", false, false);
		root.setResourceBase("/tmp/test-root");
		root.addServlet(new ServletHolder(new UploadServlet()), "/*");
		jettyInstance.start();
	}
	
	public static class UploadServlet extends DefaultServlet {
		
		@Override
		protected void doGet(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {
			System.out.println("GET request for " + request.getRequestURL());
			
			super.doGet(request, response);	
		}
		
		@Override
		protected void doPut(HttpServletRequest request, HttpServletResponse response) throws ServletException, IOException {

			System.out.println("PUT request for " + request.getRequestURL());

			
			File file = new File(getServletContext().getRealPath(URIUtil.addPaths(request.getServletPath(), request.getPathInfo())));
			if (file.exists()) {			
				// replace file if it exists
				boolean success = file.delete(); 
				if (!success) {
					response.sendError(HttpURLConnection.HTTP_INTERNAL_ERROR); // file existed and could not be deleted
					return;
				}
			}
			
			FileOutputStream out = new FileOutputStream(file);
			try {
				IO.copy(request.getInputStream(), out);
				
			} catch (IOException e) {
				//Log.warn(Log.EXCEPTION, e); // is this obsolete?
				out.close();
				throw(e);
			}
			
			response.setStatus(HttpURLConnection.HTTP_NO_CONTENT); // we return no content
		}
	}
}
