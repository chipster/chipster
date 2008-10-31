package fi.csc.microarray.util;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Reader;
import java.io.Writer;
import java.net.HttpURLConnection;

/**
 * 
 * @author Aleksi Kallio
 *
 */
public class IOUtils {

	/**
	 * Closes Reader if it is not null. Ignores all exceptions.
	 * Useful for those finally-blocks.
	 */
	public static void closeIfPossible(Reader reader) {
		if (reader != null) {
			try {
				reader.close();
			} catch (IOException e) {
				// ignore
			}
		}		
	}

	/**
	 * Closes Writer if it is not null. Ignores all exceptions.
	 * Useful for those finally-blocks.
	 */
	public static void closeIfPossible(Writer writer) {
		if (writer != null) {
			try {
				writer.close();
			} catch (IOException e) {
				// ignore
			}
		}		
	}
	
	/**
	 * Closes InputStream if it is not null. Ignores all exceptions.
	 * Useful for those finally-blocks.
	 */
	public static void closeIfPossible(InputStream stream) {
		if (stream != null) {
			try {
				stream.close();
			} catch (IOException e) {
				// ignore
			}
		}
	}
	
	/**
	 * Closes OutputStream if it is not null. Ignores all exceptions.
	 * Useful for those finally-blocks.
	 */
	public static void closeIfPossible(OutputStream stream) {
		if (stream != null) {
			try {
				stream.close();
			} catch (IOException e) {
				// ignore
			}
		}
	}

	public static void disconnectIfPossible(HttpURLConnection connection) {
		if (connection != null) {
				connection.disconnect();
		}		
	}
}
