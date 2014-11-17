package fi.csc.microarray.util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.security.KeyStore;
import java.security.KeyStoreException;
import java.security.NoSuchAlgorithmException;
import java.security.cert.CertificateException;

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
	 * @throws MicroarrayException 
	 */
	public static void initialiseTrustStore()
			throws NoSuchAlgorithmException, CertificateException,
			FileNotFoundException, IOException, KeyStoreException {
		
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();

		String password = configuration.getString("security", "storepass");
		String trustStore;
		
		if (DirectoryLayout.getInstance().getType() == DirectoryLayout.Type.CLIENT) {
			trustStore = getClientTrustStore(configuration, password);
		} else {
			//server
			trustStore = configuration.getString("security", "server-truststore");	
		}
		
		if (trustStore != null) {

			// Configure trust store
			System.setProperty("javax.net.ssl.keyStorePassword", new String(
					password));
			System.setProperty("javax.net.ssl.trustStore", trustStore);
		}
	}

	/**
	 * Return SslContextFactory with keystore initialised from Chipster
	 * configuration. To be used on the server side.
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
			String keyStorePassword)
			throws NoSuchAlgorithmException, CertificateException,
			FileNotFoundException, KeyStoreException, IOException {

		// Initialise SslContextFactory to use the keystore
		SslContextFactory sslContextFactory = new SslContextFactory(keyStore);
		sslContextFactory.setKeyStorePassword(keyStorePassword);

		return sslContextFactory;
	}
}
