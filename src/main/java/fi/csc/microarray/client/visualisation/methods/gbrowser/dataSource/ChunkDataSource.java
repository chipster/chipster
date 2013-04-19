package fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource;

import java.io.ByteArrayOutputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.util.IOUtils;

/**
 * Handler for data sources that are accessed directly, meaning that they do not
 * have indexes (like tab-separated tables). Reads data in chunks and the user
 * must parse meaningful content out from the chunks.
 * 
 * @author Petri Klemelä, Aleksi Kallio
 *
 */
public class ChunkDataSource extends DataSource {

	private RandomAccessFile raFile;

	private Long length = null;

	public ChunkDataSource(URL url, Class<? extends AreaRequestHandler> requestHandler) throws FileNotFoundException, URISyntaxException {
		super(url, requestHandler);

		if (file != null) { //Initialized by super constructor if file is local
			raFile = new RandomAccessFile(file.getPath(), "r");
		}
	}

	public ChunkDataSource(URL urlRoot, String path, Class<? extends AreaRequestHandler> requestHandler)
			throws FileNotFoundException, MalformedURLException, URISyntaxException {
		super(urlRoot, path, requestHandler);

		if (file != null) { //Initialized by super constructor if file is local
			raFile = new RandomAccessFile(file.getPath(), "r");
		}
	}
	
	public int read(long filePosition, byte[] chunk) throws IOException {
		return read(filePosition, chunk, true);
	}

	/**
	 * Method for getting a range from the file. When using url, we might not get very long ranges
	 * (like one megabyte or so).
	 * 
	 * @param filePosition
	 * @param chunk
	 * @param retry Set true to do another requests when server doesn't send as meny bytes as requested, false
	 * to try just ones.  
	 * @return
	 * @throws IOException
	 */
	public int read(long filePosition, byte[] chunk, boolean retry) throws IOException {       

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

				if (retry) {

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
				}

				return bytes;

			} catch (IOException e) {
				if(e.getMessage().contains("HTTP") && e.getMessage().contains(" 416 ")) {
					//Requested Range Not Satisfiable
					//This happens often when data files have bigger coordinates than annotations, just ignore
				} else {
					throw e;
				}
			}finally {
				IOUtils.disconnectIfPossible(connection);
			}
		}
		return -1;   
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
			ByteArrayOutputStream out = null;

			try {
				connection = (HttpURLConnection)url.openConnection();

				in = connection.getInputStream();
				out = new ByteArrayOutputStream();

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
				return out.toByteArray();
			} 
			return null;
		}
	}

	public long length() throws IOException {
		if (length == null) {
			if (raFile != null) {
				length = raFile.length();

			} else {
				HttpURLConnection connection = null;
				try {
					connection = (HttpURLConnection)url.openConnection();
					// connection.getContentLength() returns int, which is not enough
					String string = connection.getHeaderField("content-length");
					if (string == null) {
						throw new IOException("content-length unavailable for " + url);
					}
					length = Long.parseLong(connection.getHeaderField("content-length"));
				} finally {
					IOUtils.disconnectIfPossible(connection);
				}       
			} 
		}
		return length;
	}

	public RandomAccessFile getFile() {
		return raFile;
	}

	public void close() {
		if (raFile != null) {
			try {
				raFile.close();
			} catch (IOException e) {
				//No problem
			}
			raFile = null;
		}
	}
}
