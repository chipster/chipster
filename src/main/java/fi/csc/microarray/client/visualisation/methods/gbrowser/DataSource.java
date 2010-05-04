package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.RandomAccessFile;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;

import fi.csc.microarray.util.IOUtils;

public class DataSource {

	private RandomAccessFile file = null;
	
	//public for debug
	public URL url = null;
	
	public DataSource(URL url) throws FileNotFoundException {
		this.url = url;
	}

	public DataSource(File file) throws FileNotFoundException {
		this.file = new RandomAccessFile(file, "r");
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

				if(url.getFile().endsWith("seq.tsv")) {
				System.out.println(url);
				System.out.println(new String(chunk).substring(0, new String(chunk).indexOf('\n')));
				}
				
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
				int length = connection.getContentLength();
				return length;
			} finally {
				IOUtils.disconnectIfPossible(connection);
			}
		}
	}
}
