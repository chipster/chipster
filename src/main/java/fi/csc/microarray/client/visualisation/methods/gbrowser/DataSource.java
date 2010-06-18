package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;

import fi.csc.microarray.util.IOUtils;

/**
 * One source of genomic content, typically corresponding to one track. Abstraction hides
 * the physical data source, that can be either a file or a URL accessed via HTTP.
 *
 */
public class DataSource {

	private RandomAccessFile file = null;
	private URL url = null;
	private String name;
	
	public DataSource(URL url) throws FileNotFoundException {
		this.url = url;
		this.name = url.toString(); 
	}

	public DataSource(File file) throws FileNotFoundException {
		this.file = new RandomAccessFile(file, "r");
		this.name = file.toString();
	}
	
	public DataSource(URL urlRoot, String path) throws FileNotFoundException, MalformedURLException {
		this(new URL(urlRoot.toString() + "/" + path));
	}

	public DataSource(File fileRoot, String path) throws FileNotFoundException, MalformedURLException {
		this(new File(fileRoot, path));
	}

	public int read(long filePosition, byte[] chunk) throws IOException {		
		
		if (file != null) {
			file.seek(filePosition);
			return file.read(chunk);
			
		} else {
			
			HttpURLConnection connection = null;
			try {
				
				connection = (HttpURLConnection)url.openConnection();
				connection.setRequestProperty("Range", "bytes=" + filePosition + "-" + (filePosition + chunk.length));
				int bytes = connection.getInputStream().read(chunk);
				
				return bytes;
				
			} finally {
				IOUtils.disconnectIfPossible(connection);
			}
		}
		
	}

	public long length() throws IOException {
		if (file != null) {
			return file.length();
			
		} else {
			HttpURLConnection connection = null;
			try {
				connection = (HttpURLConnection)url.openConnection();
				return Long.parseLong(connection.getHeaderField("content-length")); // connection.getContentLength() returns int, which is not enough
			} finally {
				IOUtils.disconnectIfPossible(connection);
			}
		}
	}
	
	@Override
	public String toString() {
		return name;
	}
}
