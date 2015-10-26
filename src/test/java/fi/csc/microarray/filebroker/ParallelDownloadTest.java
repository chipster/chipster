package fi.csc.microarray.filebroker;



import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.Writer;
import java.net.URL;
import java.net.URLConnection;
import java.security.KeyManagementException;
import java.security.NoSuchAlgorithmException;
import java.security.cert.X509Certificate;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import javax.net.ssl.HostnameVerifier;
import javax.net.ssl.SSLSession;
import javax.net.ssl.X509TrustManager;

import org.eclipse.jetty.server.Connector;
import org.eclipse.jetty.server.Server;
import org.eclipse.jetty.server.ServerConnector;
import org.eclipse.jetty.servlet.DefaultServlet;
import org.eclipse.jetty.servlet.ServletContextHandler;
import org.eclipse.jetty.servlet.ServletHolder;
import org.eclipse.jetty.util.ssl.SslContextFactory;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;
import fi.csc.microarray.util.KeyAndTrustManager;

public class ParallelDownloadTest {
	
	public static int requests = 1000;
	public static int threads = 10;
	public static boolean readOnlyStart = true;
	
	public volatile static int done = 0;
//	private static SSLSocketFactory socketFactory;
			
	public static void main(String args[]) throws Exception {
		
		//System.setProperty("org.eclipse.jetty.LEVEL", "DEBUG");
		
//		File file = createFile();
//		final URL url = new URL("https://localhost:8080/" + file.getName());
		final URL url = new URL("https://vm0180.kaj.pouta.csc.fi:8080/storage/fe4e3b6e-331b-40a3-a485-41217a99e8e9");
//		Server server = getJetty();
//		server.start();
		
		System.out.println("Make " + requests + " requests");
				
		ExecutorService executor = Executors.newFixedThreadPool(threads);
		long t = System.currentTimeMillis();
		for (int i = 0; i < requests; i++) {
			executor.execute(new Runnable() {
				@Override
				public void run() {
					try {
						download(url);
					} catch (IOException | KeyManagementException | NoSuchAlgorithmException e) {
						e.printStackTrace();
					}
				}
			});
		}
		
		executor.shutdown();
		
		while (!executor.awaitTermination(3, TimeUnit.SECONDS)) {
			System.out.println(done + " requests completed");
		}
		System.out.println((System.currentTimeMillis() - t) + "ms");
		
		System.out.println("Wait 30 seconds for errors");
		Thread.sleep(30_000);
		System.out.println("Done");
		
//		server.stop();
//		file.delete();
	}

	public static void download(URL url) throws IOException, NoSuchAlgorithmException, KeyManagementException {
		URLConnection connection = url.openConnection();

		KeyAndTrustManager.configureForTrustAllCertificates(connection);
//		SSLContext ctx = SSLContext.getInstance("TLS");
//		ctx.init(null, new TrustManager[] { new TrustAllX509TrustManager() } , null);
//		((HttpsURLConnection) connection).setSSLSocketFactory(ctx.getSocketFactory());
//		((HttpsURLConnection) connection).setHostnameVerifier(new TrustAllHostnameVerifier());
		try {
			Thread.sleep(1);
		} catch (InterruptedException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// disable keep-alive
		//connection.setRequestProperty("connection","close");
		
		try (BufferedReader reader = new BufferedReader(new InputStreamReader(connection.getInputStream()))) {
			reader.readLine();
		}
		
		done++;
	}
	
	public static Server getJetty() throws Exception {
		Server server = new Server();
		
		//http
		//ServerConnector connector = new ServerConnector(server);        
		
		//https
		SslContextFactory factory = new SslContextFactory();
		factory.setKeyStorePath("filebroker.ks");
		factory.setKeyStorePassword("password");
		ServerConnector connector = new ServerConnector(server, factory);
		
		connector.setPort(8080);
		server.setConnectors(new Connector[]{ connector });
		
		ServletContextHandler handler = new ServletContextHandler(server, "/", false, false);
		handler.setResourceBase(new File("").getAbsolutePath());
		handler.addServlet(new ServletHolder(new DefaultServlet()), "/*");
		
		return server;
	}
	
	public static File createFile() throws IOException {
		File file = File.createTempFile("jetty-test", ".tmp", new File("").getAbsoluteFile());
		try (Writer writer = new BufferedWriter(new FileWriter(file))) {
			for (int i = 0; i < 100_000; i++) {
				writer.write("xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n");
			}
		}
		return file;
	}
	
	public static class TrustAllX509TrustManager implements X509TrustManager {
	    public X509Certificate[] getAcceptedIssuers() {
	        return new X509Certificate[0];
	    }

	    public void checkClientTrusted(java.security.cert.X509Certificate[] certs,
	            String authType) {
	    	throw new NotImplementedException();
	    }

	    public void checkServerTrusted(java.security.cert.X509Certificate[] certs,
	            String authType) {
	    }
	}
	
	public static class TrustAllHostnameVerifier implements HostnameVerifier {
	    public boolean verify(String string,SSLSession ssls) {
	        return true;
	    }
	}
}
