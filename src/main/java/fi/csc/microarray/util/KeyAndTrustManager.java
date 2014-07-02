package fi.csc.microarray.util;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.security.KeyStore;
import java.security.KeyStoreException;
import java.security.NoSuchAlgorithmException;
import java.security.cert.CertificateException;

import org.apache.log4j.Logger;
import org.eclipse.jetty.util.ssl.SslContextFactory;

public class KeyAndTrustManager {

	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(KeyAndTrustManager.class);
	
	private static void initialiseKeystore(String keyStore, String keyStorePassword, String keyAlias, String originalKeyStore)  throws NoSuchAlgorithmException, CertificateException, FileNotFoundException, IOException, KeyStoreException {
		
		// Check if proper keystore file exists
		if (!new File(keyStore).exists()) {
			logger.info("keystore file missing, exporting it");
			KeyStore original = KeyStore.getInstance(KeyStore.getDefaultType());
			original.load(KeyAndTrustManager.class.getResourceAsStream(originalKeyStore), keyStorePassword.toCharArray());
			FileOutputStream out = new FileOutputStream(keyStore); 
			original.store(out, keyStorePassword.toCharArray());
			out.close();
		}
		
		// Check if proper key exist in store
		KeyStore ks = KeyStore.getInstance(KeyStore.getDefaultType());
		FileInputStream in = new FileInputStream(keyStore);
		try {
			ks.load(in, keyStorePassword.toCharArray());
			if (!ks.containsAlias(keyAlias)) {
				logger.error("keystore is broken, bailing out");
				throw new KeyStoreException("keystore does not contain alias " + keyAlias);
			}
		} finally {
			in.close();
		}
	}

	/**
	 * Configure system to use SSL. Sets system wide properties so that Chipster key store is used. Initialises key store if needed.
	 */
	public static void initialiseSystem(String keyStore, String keyStorePassword, String keyAlias, String originalKeyStore) throws NoSuchAlgorithmException, CertificateException, FileNotFoundException, IOException, KeyStoreException {

		// Check key store and initialise if needed
		initialiseKeystore(keyStore, keyStorePassword, keyAlias, originalKeyStore);
		
		// Configure key store
		System.setProperty("javax.net.ssl.keyStore", keyStore);
		System.setProperty("javax.net.ssl.keyStorePassword", new String(keyStorePassword));
		System.setProperty("javax.net.ssl.trustStore", keyStore);						
	}


	/**
	 * Return SslContextFactory with keystore initialised from Chipster configuration. 
	 * SslContextFactory is set to trust all keys. Initialises key store if needed.
	 * 
	 * @return SslContextFactory instance
	 * @throws IOException 
	 * @throws KeyStoreException 
	 * @throws FileNotFoundException 
	 * @throws CertificateException 
	 * @throws NoSuchAlgorithmException 
	 */
	public static SslContextFactory createSslContextFactory(String keyStore, String keyStorePassword, String keyAlias, String originalKeyStore) throws NoSuchAlgorithmException, CertificateException, FileNotFoundException, KeyStoreException, IOException {
		// Check key store and initialise if needed
		initialiseKeystore(keyStore, keyStorePassword, keyAlias, originalKeyStore);

		// Initialise SslContextFactory to use the keystore
		SslContextFactory sslContextFactory = new SslContextFactory(keyStore);
		sslContextFactory.setKeyStorePassword(keyStorePassword);

		// Make SslContextFactory trust all keys
		sslContextFactory.setTrustAll(true);
		
		return sslContextFactory;
	}
}
