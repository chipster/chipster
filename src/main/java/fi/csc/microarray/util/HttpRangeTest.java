package fi.csc.microarray.util;
import java.io.ByteArrayOutputStream;
import java.io.InputStream;
import java.net.HttpURLConnection;
import java.net.URL;

import org.eclipse.jetty.server.Connector;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.ServerConnector;
import org.eclipse.jetty.servlet.DefaultServlet;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;
import org.eclipse.jetty.util.ssl.SslContextFactory;

/**
 * 
 * Setup jetty and download a part of the file
 * 
 * Doesn't really test anything yet, but provides a playground for testing 
 * Jetty and HttpURLConnection configuration.
 * 
 * @author klemela
 *
 */
public class HttpRangeTest {
	
	private static String protocol = "https";
	
	public static void main(String args[]) throws Exception {
		try {
			jetty();
			get();
		} catch (Exception e) {
			e.printStackTrace();
		} finally {
			System.exit(0);
		}
	}

	private static void jetty() throws Exception {
		System.setProperty("org.eclipse.jetty.LEVEL", "DEBUG");
		
		Server jettyInstance = new Server();
		
		ServerConnector connector = null;
		switch (protocol) {
		case "http":
			connector = new ServerConnector(jettyInstance);        
			break;
			
		case "https":			
			SslContextFactory contextFactory = KeyAndTrustManager.createSslContextFactory("../chipster-environment/security/filebroker.ks", "password", new String[] {"TLSv1.2"});
			connector = new ServerConnector(jettyInstance, contextFactory);
			break;		
		}

		connector.setPort(8080);
		jettyInstance.setConnectors(new Connector[]{ connector });

		ServletContextHandler root = new ServletContextHandler(jettyInstance, "/", false, false);
		root.setResourceBase(".");
		root.addServlet(new ServletHolder(new DefaultServlet()), "/*");
		jettyInstance.start();	
	}

	private static void get() throws Exception {
	    //for localhost testing only
	    javax.net.ssl.HttpsURLConnection.setDefaultHostnameVerifier(
	    new javax.net.ssl.HostnameVerifier(){
 
	        public boolean verify(String hostname,
	                javax.net.ssl.SSLSession sslSession) {
	            if (hostname.equals("localhost")) {
	                return true;
	            }
	            return false;
	        }
	    });
		
		String trustStore = "../chipster-environment/security/client.ts";
		String password = "password";
		
		System.setProperty("javax.net.ssl.keyStorePassword", password);
		System.setProperty("javax.net.ssl.trustStore", trustStore);
		
		URL url = new URL(protocol + "://localhost:8080/src/test/resources/affy_example.cel");
		
		long start = System.currentTimeMillis();
		
		HttpURLConnection connection = (HttpURLConnection)url.openConnection();
		connection.setRequestProperty("Range", "bytes=100000-200000");
		InputStream stream = connection.getInputStream();
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		
		IOUtils.copy(stream, out);
		
		System.out.print(new String(out.toByteArray()));
				
		long end = System.currentTimeMillis();
		
		System.out.println();
		System.out.println((end - start) + "ms");		
		System.out.println(out.size() + " bytes");		
	}
}

