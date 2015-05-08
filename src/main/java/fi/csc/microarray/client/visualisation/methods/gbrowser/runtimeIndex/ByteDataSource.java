package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.RandomAccessFile;
import java.net.HttpURLConnection;
import java.net.URISyntaxException;
import java.nio.channels.Channels;
import java.nio.channels.FileChannel;

import org.eclipse.jetty.util.IO;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.util.IOUtils;
import fi.csc.microarray.util.KeyAndTrustManager;

/**
 * Handler for data sources that are accessed directly, meaning that they do not
 * have indexes (like tab-separated tables). Reads data to byte array and the user
 * must parse meaningful content out from the bytes.
 * 
 * @author Petri Klemelä, Aleksi Kallio
 *
 */
public class ByteDataSource extends DataSource {

	private FileChannel fileChannel;
	RandomAccessFile raFile;

	private Long length = null;

	public ByteDataSource(DataUrl dataUrl) throws URISyntaxException, IOException {
		super(dataUrl);
		
		if (file != null) { //Initialized by super constructor if file is local
			raFile = new RandomAccessFile(file.getPath(), "r");
			fileChannel = raFile.getChannel();
		}
	}
	
	/**
	 * Method for getting a range from the file.
	 * 
	 * @param filePosition
	 * @param length
	 * @return array of bytes
	 * @throws IOException
	 */
	public byte[] read(long filePosition, long length) throws IOException {

		if (fileChannel != null) {
			InputStream in = Channels.newInputStream(fileChannel.position(filePosition));
			try (ByteArrayOutputStream out = new ByteArrayOutputStream()) {
				
				IO.copy(in, out, length);

				return out.toByteArray();
			}

		} else {

			long endFilePosition = filePosition + length - 1;

			//Make sure that we won't make requests outside the file end   	        	        	
			if (endFilePosition > length()) {
				endFilePosition = length();
			}

			HttpURLConnection connection = null;
			try {

				connection = (HttpURLConnection)url.openConnection();
				KeyAndTrustManager.configureForChipsterCertificate(connection);
				connection.setRequestProperty("Range", "bytes=" + filePosition + "-" + endFilePosition);
				
				try (InputStream in = connection.getInputStream();				
				ByteArrayOutputStream out = new ByteArrayOutputStream()) {
					
					IOUtils.copy(in, out);
					return out.toByteArray();
				}
				
			} catch (IOException e) {
				if(e.getMessage().contains("HTTP") && e.getMessage().contains(" 416 ")) {
					//Requested Range Not Satisfiable
					//This happens often when data files have bigger coordinates than annotations, just ignore
				} else {
					throw e;
				}
			} finally {
				IOUtils.disconnectIfPossible(connection);
			}
		}
		return null;   
	}

	/**
	 * Get all bytes from the file. Obviously this shouldn't be used for huge files, because
	 * all the data is read to the RAM. 
	 * 
	 * @return
	 */
	public byte[] readAll() throws IOException {

		return read(0, length());
	}

	public long length() throws IOException {
		if (length == null) {
			if (fileChannel != null) {
				length = file.length();

			} else {
				HttpURLConnection connection = null;
				try {
					connection = (HttpURLConnection)url.openConnection();
					KeyAndTrustManager.configureForChipsterCertificate(connection);
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

	public void close() {
		if (raFile != null) {
			try {
				raFile.close();
			} catch (IOException e) {
				//No problem
			}
			fileChannel = null;
			raFile = null;
		}
	}
}
