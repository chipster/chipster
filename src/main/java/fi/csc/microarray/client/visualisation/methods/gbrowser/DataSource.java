package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.RandomAccessFile;
import java.net.HttpURLConnection;
import java.net.MalformedURLException;
import java.net.URL;

import fi.csc.microarray.util.IOUtils;

public class DataSource {

	private RandomAccessFile file = null;
	private HttpURLConnection connection = null;
	private BufferedReader urlReader = null;
	private long urlReaderPosition = -1; 
	
	public DataSource(URL url) throws FileNotFoundException, IOException {
		
		connection = (HttpURLConnection)url.openConnection();
	}

	public DataSource(File file) throws FileNotFoundException {
		this.file = new RandomAccessFile(file, "r");
	}
	
	public DataSource(URL urlRoot, String path) throws FileNotFoundException, MalformedURLException, IOException {
		this(new URL(urlRoot.toString() + "/" + path));
	}

	public DataSource(File fileRoot, String path) throws FileNotFoundException, MalformedURLException {
		this(new File(fileRoot, path));
	}
	
	public void seek(long filePosition) throws IOException {
		if (file != null) {
			file.seek(filePosition);
			
		} else {
			long length = connection.getContentLength();
			
			connection.disconnect();
			connection.setRequestProperty("Range", "bytes=" + filePosition + "-" + length);
			urlReader = new BufferedReader(new InputStreamReader(connection.getInputStream()));
			urlReaderPosition = filePosition;
		}		
	}
	
	public InputStream getStream() throws IOException {
		if (file != null) {
			throw new IllegalStateException("Can't get input stream from random access file");
			
		} else {
			//connection.setRequestProperty("Range", "bytes=" + filePosition + "-" + connection.getContentLength());
			return connection.getInputStream();
		}		
	}

	public int read(byte[] chunk) throws IOException {
		
		if (file != null) {
			return file.read(chunk);
			
		} else {
			if (connection == null) {
				throw new IllegalStateException("Seek has to be called first");
			}
			int bytes = connection.getInputStream().read(chunk);
			urlReaderPosition += bytes;
			return bytes;
		}		
	}
	
	public String readLine() throws IOException {
		
		if (file != null) {
			return file.readLine();
			
		} else {
			
			String line =  urlReader.readLine();
			
			//plus one to take care of the new line character removed by readLine()
			if (line != null) {
				urlReaderPosition += line.length();
				
				if (line.charAt(line.length() - 1) == '\n') {
					urlReaderPosition++;
				}
			}
			
			return line;
		}		
	}
	
	public long getPosition() throws IOException {
		
		if (file != null) {
			return file.getFilePointer();
		} else {
			return urlReaderPosition;
		}
	}
	
	public void clean() throws IOException {
		if (file != null) {

		} else {
			urlReader.close();
			IOUtils.disconnectIfPossible(connection);
		}		
	}

	public long length() throws IOException {
		if (file != null) {
			return file.length();
			
		} else {
			int length = connection.getContentLength();
			return length;			
		}
	}
}
