package fi.csc.microarray.util;

import java.io.BufferedReader;
import java.io.IOException;

/**
 * 
 * @author Aleksi Kallio
 *
 */
public class IOUtils {

	/**
	 * Closes reader if it is not null. Ignores all exceptions.
	 * Useful for those finally-blocks.
	 */
	public static void closeIfPossible(BufferedReader reader) {
		if (reader != null) {
			try {
				reader.close();
			} catch (IOException e) {
				// ignore
			}
		}
		
	}

}
