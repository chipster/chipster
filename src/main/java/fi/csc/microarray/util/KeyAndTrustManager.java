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

public class KeyAndTrustManager {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(KeyAndTrustManager.class);
	
	
	public static void initialise(String keyStore, char[] keyStorePassword, String keyAlias, String originalKeyStore) throws NoSuchAlgorithmException, CertificateException, FileNotFoundException, IOException, KeyStoreException {
		
		// check if proper keystore file exists
		if (!new File(keyStore).exists()) {
			logger.info("keystore file missing, exporting it");
			KeyStore original = KeyStore.getInstance(KeyStore.getDefaultType());
			original.load(KeyAndTrustManager.class.getResourceAsStream(originalKeyStore), keyStorePassword);
			FileOutputStream out = new FileOutputStream(keyStore); 
			original.store(out, keyStorePassword);
			out.close();
		}
		
		// check if proper keys exist
		KeyStore ks = KeyStore.getInstance(KeyStore.getDefaultType());
		FileInputStream in = new FileInputStream(keyStore);
		try {
			ks.load(in, keyStorePassword);
			if (!ks.containsAlias(keyAlias)) {
				logger.error("keystore is broken, bailing out");
				throw new KeyStoreException("keystore does not contain alias " + keyAlias);
			}
		} finally {
			in.close();
		}

		System.setProperty("javax.net.ssl.keyStore", keyStore);
		System.setProperty("javax.net.ssl.keyStorePassword", new String(keyStorePassword));
		System.setProperty("javax.net.ssl.trustStore", keyStore);
	}
}
