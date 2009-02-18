package fi.csc.microarray.security;

import java.util.UUID;

/**
 * Cryptographically strong pseudorandom key. Uses Java UUID for actual implementation.
 *  
 * @author Aleksi Kallio
 *
 */
public class CryptoKey {

	public static String generateRandom() {
		return UUID.randomUUID().toString();
	}
	
	public static boolean validateKeySyntax(String string) {
		try {
			UUID.fromString(string);
			return true; // parse successful
			
		} catch (IllegalArgumentException e) {
			return false; // parse failed
		}
	}
}
