package fi.csc.microarray.util;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
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

	private static final int BUFFER_SIZE = 16*1024;
	private static final long CALLBACK_INTERVAL = 500; 
	
	/**
	 * Closes Reader if it is not null. Ignores all exceptions. Useful for those finally-blocks.
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
	 * Closes Writer if it is not null. Ignores all exceptions. Useful for those finally-blocks.
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
	 * Closes InputStream if it is not null. Ignores all exceptions. Useful for those finally-blocks.
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
	 * Closes OutputStream if it is not null. Ignores all exceptions. Useful for those finally-blocks.
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

	public static interface CopyProgressListener {
		public void progress(int bytes);
	}
	
	/**
	 * 
	 * Copies stream contents of source to target and reports progress.
	 * 
	 * @param source input stream
	 * @param target output stream
	 * @param progressListener can be null
	 * 
	 * @throws IOException all exceptions from underlying IO are passed through
	 */
	public static void copy(InputStream source, OutputStream target, CopyProgressListener progressListener) throws IOException {
		
		// initialise
		byte buffer[] = new byte[BUFFER_SIZE];
		int len = BUFFER_SIZE;
		int sum = 0;
		long lastCallback = Long.MAX_VALUE; 
		
		// tell that we are in the beginning
		if (progressListener != null) {
			progressListener.progress(0);
			lastCallback = System.currentTimeMillis();
		}

		// copy while there is content
		while (true) {
			len = source.read(buffer, 0, BUFFER_SIZE);
			if (len < 0) {
				break;
			}
			
			target.write(buffer, 0, len);
			sum += len;
			
			// report progress every CALLBACK_INTERVAL milliseconds
			if (progressListener != null && (lastCallback+CALLBACK_INTERVAL) < System.currentTimeMillis()) {
				progressListener.progress(sum);
				lastCallback = System.currentTimeMillis();
			}
		}
	}
	
	public static void copy(InputStream source, OutputStream target) throws IOException {
		copy(source, target, null);
	}

	/**
	 * Copies a file. 
	 * 
	 * @param from source file
	 * @param to destination file
	 * 
	 * @throws IOException if copying if file contents fails
	 */
	public static void copy(File from, File to) throws IOException {
		FileInputStream in = new FileInputStream(from);
		FileOutputStream out = new FileOutputStream(to);
		
		try {
			copy(in, out, null);
			
		} finally {
			closeIfPossible(in);
			closeIfPossible(out);
		}
	}

}
