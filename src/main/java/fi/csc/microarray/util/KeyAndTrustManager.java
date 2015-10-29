package fi.csc.microarray.util;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.net.URL;
import java.net.URLConnection;
import java.security.KeyManagementException;
import java.security.KeyStore;
import java.security.KeyStoreException;
import java.security.NoSuchAlgorithmException;
import java.security.cert.CertificateException;
import java.security.cert.X509Certificate;

import javax.net.ssl.HostnameVerifier;
import javax.net.ssl.HttpsURLConnection;
import javax.net.ssl.SSLContext;
import javax.net.ssl.SSLSession;
import javax.net.ssl.SSLSocketFactory;
import javax.net.ssl.TrustManager;
import javax.net.ssl.TrustManagerFactory;
import javax.net.ssl.X509TrustManager;

import org.apache.commons.io.FileUtils;
import org.apache.log4j.Logger;
import org.eclipse.jetty.util.ssl.SslContextFactory;

import sun.reflect.generics.reflectiveObjects.NotImplementedException;
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

	/**
	 * SSLSocketFactory with a truststore containing Chipster server's
	 * self-signed certificate, or null when the server has a certificate signed
	 * by a CA.
	 */
	private static SocketFactoryCache selfSignedFactoryCache;

	/**
	 * Default SSLSocketFactory provided by Java with all the default CA root
	 * certificates.
	 */
	private static SSLSocketFactory caFactory;

	private static SocketFactoryCache trustAllFactoryCache;

	private static String protocol = "TLS";

	public static String getClientTrustStore(Configuration configuration,
			String password) throws NoSuchAlgorithmException,
			CertificateException, FileNotFoundException, IOException,
			KeyStoreException {

		String trustStoreFilename = configuration.getString("security", "client-truststore");

		if (trustStoreFilename == null || "".equals(trustStoreFilename)) {
			// CA signed certificate
			return null;
		}

		String brokerHost = configuration.getString("messaging", "broker-host");
		// it is assumed that a trust store is located next to the configuration
		// file on a server
		URL remoteTrustStore = new URL(configuration.getConfigRootURL()
				+ trustStoreFilename);
		// there must be a separate local file for each installation
		File localTrustStore = new File(DirectoryLayout.getInstance()
				.getSecurityDir()
				+ File.separator
				+ brokerHost
				+ "-"
				+ trustStoreFilename);

		// Check if proper keystore file exists
		if (!localTrustStore.exists()) {

			logger.info("keystore file missing, downloading it");
			// rename file only when it's complete to avoid problems when there
			// are parallel instances of SessionReplayTest
			File downloadTemp = File.createTempFile(localTrustStore.getName(), "", localTrustStore.getParentFile());
			FileUtils.copyURLToFile(remoteTrustStore, downloadTemp);
			if (!downloadTemp.renameTo(localTrustStore)) {
				logger.info("renaming keystore file failed, maybe other client is running at the same time? (error ignored)");
			}
		}

		return localTrustStore.getPath();
	}

	/**
	 * Configure system to use SSL. Sets system wide properties so that Chipster
	 * key store is used. Initialises key store if needed.
	 * 
	 * @throws KeyManagementException
	 * @throws MicroarrayException
	 */
	public static void initialiseTrustStore() throws NoSuchAlgorithmException,
			CertificateException, FileNotFoundException, IOException,
			KeyStoreException, KeyManagementException {

		if (!initialised) {

			// take a copy of the default factory accepting CA certificates
			caFactory = getSocketFactory(null);

			Configuration configuration = DirectoryLayout.getInstance().getConfiguration();

			String password = configuration.getString("security", "storepass");
			String trustStorePath;

			if (DirectoryLayout.getInstance().getType() == DirectoryLayout.Type.CLIENT) {
				trustStorePath = getClientTrustStore(configuration, password);

			} else {
				// server
				trustStorePath = configuration.getString("security", "server-truststore");
				if ("".equals(trustStorePath)) {
					trustStorePath = null;
				} else if (!new File(trustStorePath).exists()) {
					throw new RuntimeException("configured trust store file " + trustStorePath + " doesn't exist");
				}
			}

			if (trustStorePath != null) {

				// configure trust store
				// this should be enough, but see comment about sslFactory below
				System.setProperty("javax.net.ssl.keyStorePassword", password);
				System.setProperty("javax.net.ssl.trustStore", trustStorePath);

				if (DirectoryLayout.getInstance().getType() == DirectoryLayout.Type.CLIENT) {
					// init sslFactory (hides complaints about untrusted
					// certificates
					// when started with web-start)
					KeyStore trustStore = KeyStore.getInstance(KeyStore.getDefaultType());

					try (FileInputStream inputStream = new FileInputStream(trustStorePath)) {
						trustStore.load(inputStream, password.toCharArray());
					}

					TrustManagerFactory tmf = TrustManagerFactory.getInstance(TrustManagerFactory.getDefaultAlgorithm());
					tmf.init(trustStore);
					selfSignedFactoryCache = new SocketFactoryCache(tmf.getTrustManagers(), protocol);
				}
			}

			trustAllFactoryCache = new SocketFactoryCache(new TrustManager[] { new TrustAllX509TrustManager() }, protocol);

			if (!configuration.getBoolean("security", "verify-hostname")) {
				// Accept all hostnames (do not try to match certificate
				// hostname (CN) to observed hostname)
				HttpsURLConnection.setDefaultHostnameVerifier(new TrustAllHostnameVerifier());
			}

			initialised = true;
		}
	}

	private static SSLSocketFactory getSocketFactory(
			TrustManager[] trustManagers) throws NoSuchAlgorithmException,
			KeyManagementException {
		
		SSLContext ctx = SSLContext.getInstance("TLS");
		ctx.init(null, trustManagers, null);
		return ctx.getSocketFactory();
	}

	/**
	 * Return SslContextFactory with keystore initialised from Chipster
	 * configuration. To be used on the server side.
	 * 
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
	public static void configureForChipsterCertificate(URLConnection connection) {
		// only in client config
		if (selfSignedFactoryCache != null) {
			if (connection instanceof HttpsURLConnection) {
				try {
					((HttpsURLConnection) connection).setSSLSocketFactory(selfSignedFactoryCache.getSocketFactoryForThisThread());
				} catch (KeyManagementException | NoSuchAlgorithmException e) {
					// client can't fix these
					new RuntimeException(e);
				}
			}
		}
	}

	public static void configureForCACertificates(URLConnection connection) {
		if (connection instanceof HttpsURLConnection) {
			/*
			 * Sharing caFactory between many threads may cause problems (see
			 * SocketFactoryCache), but I'm not sure how to make new instances
			 * safely.
			 */
			((HttpsURLConnection) connection).setSSLSocketFactory(caFactory);
		}
	}

	public static void configureForTrustAllCertificates(URLConnection connection)
			throws NoSuchAlgorithmException, KeyManagementException {
		
		if (connection instanceof HttpsURLConnection) {
			((HttpsURLConnection) connection).setSSLSocketFactory(trustAllFactoryCache.getSocketFactoryForThisThread());
			((HttpsURLConnection) connection).setHostnameVerifier(new TrustAllHostnameVerifier());
		}
	}

	public static class TrustAllX509TrustManager implements X509TrustManager {
		public X509Certificate[] getAcceptedIssuers() {
			return new X509Certificate[0];
		}

		public void checkClientTrusted(
				java.security.cert.X509Certificate[] certs, String authType) {
			throw new NotImplementedException();
		}

		public void checkServerTrusted(
				java.security.cert.X509Certificate[] certs, String authType) {
		}
	}

	public static class TrustAllHostnameVerifier implements HostnameVerifier {
		public boolean verify(String string, SSLSession ssls) {
			return true;
		}
	}
}
