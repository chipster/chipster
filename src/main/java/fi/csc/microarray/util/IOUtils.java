package fi.csc.microarray.util;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.RandomAccessFile;
import java.io.Reader;
import java.io.Writer;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;

import javax.servlet.http.HttpServletRequest;

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

	public static void disconnectIfPossible(URLConnection connection) {
		if (connection instanceof HttpURLConnection) {
			disconnectIfPossible((HttpURLConnection)connection);
		}
	}
	

	public static interface CopyProgressListener {
		public void progress(long bytes);
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
		
		BufferedInputStream bSource = new BufferedInputStream(source);
		BufferedOutputStream bTarget = new BufferedOutputStream(target);
		
		// initialise
		byte buffer[] = new byte[BUFFER_SIZE];
		int len = BUFFER_SIZE;
		long sum = 0;
		long lastCallback = Long.MAX_VALUE; 
		
		// tell that we are in the beginning
		if (progressListener != null) {
			progressListener.progress(0);
			lastCallback = System.currentTimeMillis();
		}

		// copy while there is content
		while (true) {
			len = bSource.read(buffer, 0, BUFFER_SIZE);
			if (len < 0) {
				break;
			}
			
			bTarget.write(buffer, 0, len);
			sum += len;
			
			// report progress every CALLBACK_INTERVAL milliseconds
			if (progressListener != null && (lastCallback+CALLBACK_INTERVAL) < System.currentTimeMillis()) {
				progressListener.progress(sum);
				lastCallback = System.currentTimeMillis();
			}
		}
		bTarget.flush();
	}
	
	public static void copy(InputStream source, OutputStream target) throws IOException {
		copy(source, target, null);
	}
	
	public static void copy(InputStream source, File target) throws IOException {
		FileOutputStream out = new FileOutputStream(target);
		try {
			copy(source, out, null);
		} finally {
			closeIfPossible(out);
		}
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

	public static void closeIfPossible(RandomAccessFile raf) {
		if (raf != null) {
			try {
				raf.close();
			} catch (IOException e) {
				// ignore
			}
		}		
	}

	public static URL createURL(URL url, String postfix) throws MalformedURLException {
		return new URL(url, url.getFile() + "/" + postfix);
	}
	
	public static boolean isLocalFileURL(URL url) {
		return "file".equals(url.getProtocol());
	}
	
	public static String getFilenameWithoutPath(URL url) {
		return url.getPath().substring(url.getPath().lastIndexOf('/') + 1);
	}

	public static String getFilenameWithoutPath(HttpServletRequest request) {
		String path = request.getPathInfo();
		return path.substring(path.lastIndexOf('/') + 1);
	}

	/**
	 * Compare the contents of two Streams to determine if they are equal or not.
	 * 
	 * This method buffers the input internally using <code>BufferedInputStream</code> if they are
	 * not already buffered.
	 * 
	 * @param input1
	 *            the first stream
	 * @param input2
	 *            the second stream
	 * @return true if the content of the streams are equal or they both don't exist, false
	 *         otherwise
	 * @throws NullPointerException
	 *             if either input is null
	 * @throws IOException
	 *             if an I/O error occurs
	 */
	public static boolean contentEquals(InputStream input1, InputStream input2) throws IOException {
		if (!(input1 instanceof BufferedInputStream)) {
			input1 = new BufferedInputStream(input1);
		}
		if (!(input2 instanceof BufferedInputStream)) {
			input2 = new BufferedInputStream(input2);
		}

		int ch = input1.read();
		while (-1 != ch) {
			int ch2 = input2.read();
			if (ch != ch2) {
				return false;
			}
			ch = input1.read();
		}

		int ch2 = input2.read();
		return (ch2 == -1);
	}

	
	/**
	 * Gets filename of various object types, currently supporting: File, URL. If type is unknown, null
	 * is returned.
	 *  
	 * @param input some object with "filename" type of property
	 * @return
	 */
	public static String getFilename(Object o) {
		
		if (o instanceof File) {
			return ((File)o).getName();
		
		} else if (o instanceof URL) {
			return getFilenameWithoutPath((URL)o);

		} else {
			return null;
		}
	}

	/**
	 * Disable cache for the URLConnection class.
	 * 
	 * When using web start, http caching is enabled. This is a problem with big data files.
	 * 
	 */
	public static void disableHttpCache() {
		try {
			// no static way to do this, see http://bugs.sun.com/bugdatabase/view_bug.do?bug_id=4851466
			// openConnection does not really open the connection
			new URL("http://chipster.csc.fi").openConnection().setDefaultUseCaches(false);
		} catch (Exception e) {
			
		}
	}
}
