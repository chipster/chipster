package fi.csc.microarray.filebroker.encryption;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;

import javax.servlet.ServletException;
import javax.servlet.http.HttpServletRequest;
import javax.servlet.http.HttpServletResponse;

import org.mortbay.jetty.Connector;
import org.mortbay.jetty.Server;
import org.mortbay.jetty.nio.SelectChannelConnector;
import org.mortbay.jetty.security.SslSocketConnector;
import org.mortbay.jetty.servlet.Context;
import org.mortbay.jetty.servlet.DefaultServlet;
import org.mortbay.jetty.servlet.ServletHolder;
import org.mortbay.log.Log;
import org.mortbay.thread.QueuedThreadPool;
import org.mortbay.util.IO;
import org.mortbay.util.URIUtil;

import sun.net.www.protocol.http.HttpURLConnection;

public class TestServer {

	
	public static void main(String[] args) throws Exception {
		Server jettyInstance = new Server();
		jettyInstance.setThreadPool(new QueuedThreadPool());
		SelectChannelConnector connector = new SelectChannelConnector();
		connector.setServer(jettyInstance);
		connector.setPort(8080);
		SslSocketConnector sslConnector = new SslSocketConnector();
		sslConnector.setServer(jettyInstance);
		sslConnector.setPort(8443);
		sslConnector.setKeystore("keystore.ks");
		sslConnector.setKeyPassword("microarray");
		jettyInstance.setConnectors(new Connector[]{ connector, sslConnector });

		Context root = new Context(jettyInstance, "/", false, false);
		root.setResourceBase("/tmp");
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
				Log.warn(Log.EXCEPTION, e); // is this obsolete?
				out.close();
				throw(e);
			}
			
			response.setStatus(HttpURLConnection.HTTP_NO_CONTENT); // we return no content
		}
	}
}
