package fi.csc.microarray.util;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.net.URLConnection;
import java.security.KeyManagementException;
import java.security.KeyStore;
import java.security.KeyStoreException;
import java.security.NoSuchAlgorithmException;
import java.security.cert.CertificateException;

import javax.net.ssl.HostnameVerifier;
import javax.net.ssl.HttpsURLConnection;
import javax.net.ssl.SSLContext;
import javax.net.ssl.SSLSession;
import javax.net.ssl.SSLSocketFactory;
import javax.net.ssl.TrustManagerFactory;

import org.apache.log4j.Logger;
import org.eclipse.jetty.util.ssl.SslContextFactory;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.exception.MicroarrayException;

public class KeyAndTrustManager {

	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger
			.getLogger(KeyAndTrustManager.class);
	
	private static boolean initialised;

	private static SSLSocketFactory sslFactory;

	public static String getClientTrustStore(Configuration configuration, String password)
			throws NoSuchAlgorithmException, CertificateException,
			FileNotFoundException, IOException, KeyStoreException {

		String trustStoreFilename = configuration.getString("security", "client-truststore");
		
		if (trustStoreFilename == null || "".equals(trustStoreFilename)) {
			// CA signed certificate
			return null;
		}
		
		String brokerHost = configuration.getString("messaging", "broker-host");
		// it is assumed that a trust store is located next to the configuration file on a server
		URL remoteTrustStore = new URL(configuration.getConfigRootURL() + trustStoreFilename);
		// there must be a separate local file for each installation
		File localTrustStore = new File(DirectoryLayout.getInstance().getSecurityDir() + File.separator + brokerHost + "-" + trustStoreFilename);
		
		// Check if proper keystore file exists
		if (!localTrustStore.exists()) {
			
			// in practice just a complicated way of doing a http download
			
			logger.info("keystore file missing, exporting it");
			KeyStore original = KeyStore.getInstance(KeyStore.getDefaultType());

			try (InputStream urlStream = remoteTrustStore.openStream()){
				original.load(urlStream, password.toCharArray());
			} catch (FileNotFoundException e) {
				// there is no "cause" in FileNotFoundException
				throw new KeyStoreException("cannot load server's SSL certificate", e);
			}
			
			try (FileOutputStream out = new FileOutputStream(localTrustStore)) {
				original.store(out, password.toCharArray());
			}
		}
		
		return localTrustStore.getPath();
	}

	/**
	 * Configure system to use SSL. Sets system wide properties so that Chipster
	 * key store is used. Initialises key store if needed.
	 * @throws KeyManagementException 
	 * @throws MicroarrayException 
	 */
	public static void initialiseTrustStore()
			throws NoSuchAlgorithmException, CertificateException,
			FileNotFoundException, IOException, KeyStoreException, KeyManagementException {

		if (!initialised) {

			Configuration configuration = DirectoryLayout.getInstance().getConfiguration();

			String password = configuration.getString("security", "storepass");
			String trustStorePath;

			if (DirectoryLayout.getInstance().getType() == DirectoryLayout.Type.CLIENT) {
				trustStorePath = getClientTrustStore(configuration, password);
			
			} else {
				//server
				trustStorePath = configuration.getString("security", "server-truststore");
				if ("".equals(trustStorePath)) {
					trustStorePath = null;
				}
			}

			if (trustStorePath != null) {					
				// configure trust store
				// this should be enough, but see comment about sslFactory below
				System.setProperty("javax.net.ssl.keyStorePassword", password);
				System.setProperty("javax.net.ssl.trustStore", trustStorePath);
				
				if (DirectoryLayout.getInstance().getType() == DirectoryLayout.Type.CLIENT) {
					// init sslFactory (hides complaints about untrusted certificates
					// when started with web-start)
					KeyStore trustStore = KeyStore.getInstance(KeyStore.getDefaultType());

					try (FileInputStream inputStream = new FileInputStream(trustStorePath)) {
						trustStore.load(inputStream, password.toCharArray());
					}

					TrustManagerFactory tmf = TrustManagerFactory.getInstance(TrustManagerFactory.getDefaultAlgorithm());
					tmf.init(trustStore);
					SSLContext ctx = SSLContext.getInstance("TLS");
					ctx.init(null, tmf.getTrustManagers(), null);
					sslFactory = ctx.getSocketFactory();
				}
			}


			if (!configuration.getBoolean("security", "verify-hostname")) {
				// Accept all hostnames (do not try to match certificate hostname (CN) to observed hostname)
				HttpsURLConnection.setDefaultHostnameVerifier(new HostnameVerifier() {
					@Override
					public boolean verify(String hostname, SSLSession session) {
						return true; 
					}
				});
			}

			initialised = true;
		}
	}

	/**
	 * Return SslContextFactory with keystore initialised from Chipster
	 * configuration. To be used on the server side.
	 * @param protocols 
	 * 
	 * @return SslContextFactory instance
	 * @throws IOException
	 * @throws KeyStoreException
	 * @throws FileNotFoundException
	 * @throws CertificateException
	 * @throws NoSuchAlgorithmException
	 * @throws MicroarrayException 
	 */
	public static SslContextFactory createSslContextFactory(String keyStore,
			String keyStorePassword, String[] protocols)
			throws NoSuchAlgorithmException, CertificateException,
			FileNotFoundException, KeyStoreException, IOException {

		// Initialise SslContextFactory to use the keystore
		SslContextFactory sslContextFactory = new SslContextFactory(keyStore);
		sslContextFactory.setKeyStorePassword(keyStorePassword);		
		sslContextFactory.setIncludeProtocols(protocols);		
				
		return sslContextFactory;
	}

	/**
	 * Make URLConnection to trust file broker's self-signed certificate also
	 * when using web-start.
	 * 
	 * @param connection
	 */
	public static void configureSSL(URLConnection connection) {
		// only in client config
		if (sslFactory != null) {
			if (connection instanceof HttpsURLConnection) {
				((HttpsURLConnection)connection).setSSLSocketFactory(sslFactory);
			}
		}
	}
}
