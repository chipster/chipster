package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;

import com.sun.xml.messaging.saaj.util.ByteOutputStream;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.TsvParser;
import fi.csc.microarray.util.IOUtils;

/**
 * Handler for files accessed directly (e.g. tab-separated files).
 *
 */
public class ChunkDataSource extends DataSource {

	private RandomAccessFile raFile;
	private TsvParser fileParser;
	private Long length = null;

	public ChunkDataSource(URL url, TsvParser fileParser) throws FileNotFoundException {
		super(url);
		this.fileParser = fileParser;
	}

	public ChunkDataSource(File file, TsvParser fileParser) throws FileNotFoundException {
		super(file);
		raFile = new RandomAccessFile(file, "r");
		this.fileParser = fileParser;
	}

	public ChunkDataSource(URL urlRoot, String path, TsvParser fileParser)
	throws FileNotFoundException, MalformedURLException {
		super(urlRoot, path);
		this.fileParser = fileParser;
	}

	public ChunkDataSource(File fileRoot, String path, TsvParser fileParser)
	throws FileNotFoundException, MalformedURLException {
		this(new File(fileRoot, path), fileParser);
	}

	/**
	 * Method for getting a range from the file. When using url, we might not get very long ranges
	 * (like one megabyte or so).
	 * 
	 * @param filePosition
	 * @param chunk
	 * @return
	 * @throws IOException
	 */
	public int read(long filePosition, byte[] chunk) throws IOException {       

		if (raFile != null) {
			raFile.seek(filePosition);
			return raFile.read(chunk);

		} else {

			long endFilePosition = filePosition + chunk.length;

			//The ChunkFileFetcherThread makes requests bigger to get the last line fully
			//Make sure that we won't make requests outside the file end   	        	        	
			if (endFilePosition > length()) {
				endFilePosition = length();
			}

			HttpURLConnection connection = null;
			try {

				connection = (HttpURLConnection)url.openConnection();
				connection.setRequestProperty("Range", "bytes=" + filePosition + "-" + endFilePosition);
				int bytes = connection.getInputStream().read(chunk);


				/* reading seems to give only the beginning of the range (e.g. 24576 bytes) if 
				 * too big range is requested. Here we try to read again the missing part of the 
				 * range, recursively if needed. If we didn't get any bytes, there is probably no
				 * point to continue.
				 */

				if (bytes < endFilePosition - filePosition && bytes != 0) {

					byte[] nextBytes = new byte[(int) (endFilePosition - filePosition - bytes)];

					int nextLength = read(filePosition + bytes, nextBytes);

					for (int i = 0; i < nextLength; i++) {
						chunk[bytes + i] = nextBytes[i];
					}
				}

				return bytes;

			} finally {
				IOUtils.disconnectIfPossible(connection);
			}
		}   
	}

	/**
	 * Get all bytes from the file. Obviously this shouldn't be used for huge files, because
	 * all the data is read to the RAM. The same data can be obtained with method read(), but
	 * it seems to have difficulties with data sizes above a megabyte. 
	 * 
	 * @return
	 */
	public byte[] readAll() throws IOException {

		if (raFile != null) {

			byte[] buf = null;

			raFile.seek(0);				
			buf = new byte[(int) raFile.length()];
			raFile.read(buf);

			return buf;

		} else {

			HttpURLConnection connection = null;    	
			InputStream in = null;
			ByteOutputStream out = null;

			try {
				connection = (HttpURLConnection)url.openConnection();

				in = connection.getInputStream();
				out = new ByteOutputStream();

				IOUtils.copy(in, out);

			} finally {

				if (in != null) {
					in.close();
				}

				if (out != null) {
					out.close();
				}

				IOUtils.disconnectIfPossible(connection);
			}

			if (out != null) {
				return out.getBytes();
			} 
			return null;
		}
	}

	public long length() throws IOException {
		if (raFile != null) {
			return raFile.length();

		} else {
			if (length == null) {
				HttpURLConnection connection = null;
				try {
					connection = (HttpURLConnection)url.openConnection();
					// connection.getContentLength() returns int, which is not enough
					length = Long.parseLong(connection.getHeaderField("content-length"));
				} finally {
					IOUtils.disconnectIfPossible(connection);
				}
			} 
			return length;
		}
	}

	public long headerLength() {
		return fileParser.getHeaderLength(file);
	}

	public TsvParser getFileParser() {
		return fileParser;
	}


	public RandomAccessFile getFile() {
		return raFile;
	}
}
