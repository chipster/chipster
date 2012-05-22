package fi.csc.microarray.databeans.handlers;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URL;

import fi.csc.microarray.databeans.DataBean.ContentLocation;

public class RemoteContentHandler implements ContentHandler {
	
	@Override
	public InputStream getInputStream(ContentLocation location) throws FileNotFoundException {
		checkCompatibility(location);
		return null; // FIXME
	}

	@Override
	public OutputStream getOutputStream(ContentLocation location) throws IOException {
		checkCompatibility(location);
		return null; // FIXME
	}
	
	@Override
	public long getContentLength(ContentLocation location) {
		checkCompatibility(location);
		
		return 0; // FIXME
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

}
