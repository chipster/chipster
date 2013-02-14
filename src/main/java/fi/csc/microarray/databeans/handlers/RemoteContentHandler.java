package fi.csc.microarray.databeans.handlers;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;

import fi.csc.microarray.databeans.DataBean.ContentLocation;
import fi.csc.microarray.util.IOUtils;

public class RemoteContentHandler implements ContentHandler {
	
	@Override
	public InputStream getInputStream(ContentLocation location) throws IOException {
		checkCompatibility(location);
		HttpURLConnection connection = (HttpURLConnection)location.getUrl().openConnection();
		return connection.getInputStream();
	}

	@Override
	public OutputStream getOutputStream(ContentLocation location) throws IOException {
		throw new UnsupportedOperationException("remote content handler does not support output");
	}
	
	@Override
	public long getContentLength(ContentLocation location) throws IOException {
		checkCompatibility(location);
		
		HttpURLConnection connection = null;
		try {
			connection = (HttpURLConnection)location.getUrl().openConnection();
			return Long.parseLong(connection.getHeaderField("content-length"));
		} finally {
			IOUtils.disconnectIfPossible(connection);
		}
	}

	/**
	 * Deletes files that have been created for this bean. Only deletes temporary local files,
	 * not user files etc.
	 */
	@Override
	public void markDeletable(ContentLocation location) {
		// remote data is not currently deleted
	}
	
	@Override
	public void checkCompatibility(ContentLocation location) throws IllegalArgumentException {
		
		URL url = location.getUrl();
		
		// null url
		if (url == null) {
			throw new IllegalArgumentException("url is null.");
		} 

		// protocol not HTTP compatible
		else if (!("http".equals(url.getProtocol()) || "https".equals(url.getProtocol()))) {
			throw new IllegalArgumentException("Protocol of " + url.toString() + " is not HTTP compatible.");
		}
		
		// null or empty path
		else if (url.getPath() == null || url.getPath().length() == 0) {
			throw new IllegalArgumentException("Illegal path: " + url.toString());
		} 


	}

	@Override
	public boolean isAccessible(ContentLocation location) {
		checkCompatibility(location);
		try {
			HttpURLConnection connection = (HttpURLConnection)location.getUrl().openConnection();
			connection.connect () ; 
			return connection.getResponseCode() == 200;
		} catch (IOException e) {
			return false;
		}
	}

}
