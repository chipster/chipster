package fi.csc.microarray.databeans.handlers;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.HashMap;
import java.util.Map;

import org.apache.log4j.Logger;

import de.schlichtherle.truezip.zip.ZipEntry;
import de.schlichtherle.truezip.zip.ZipFile;
import fi.csc.microarray.client.session.UserSession;
import fi.csc.microarray.databeans.DataBean.ContentLocation;

public class ZipContentHandler implements ContentHandler {
	/**
	 * Logger for this class
	 */
	private static final Logger logger = Logger.getLogger(ZipContentHandler.class);

	private Map<File, ZipFile> zipFileInstances = new HashMap<File, ZipFile>();
	
	@Override
	public long getContentLength(ContentLocation location) throws IOException {
		checkCompatibility(location);
		ZipFile zipFile = createZipFile(location);
		ZipEntry zipEntry = zipFile.getEntry(location.getUrl().getRef());
		return zipEntry.getSize();
	}

	@Override
	public InputStream getInputStream(ContentLocation location) throws IOException {
		checkCompatibility(location);
		ZipFile zipFile = createZipFile(location);
		ZipEntry zipEntry = zipFile.getEntry(location.getUrl().getRef());
		return zipFile.getInputStream(zipEntry);
	}

	/**
	 * Pools ZipFile instances to avoid having too many files open.
	 * 
	 */
	private ZipFile createZipFile(ContentLocation location) throws IOException {
		File file = getZipFile(location);
		if (!zipFileInstances.containsKey(file)) {
			zipFileInstances.put(file, new ZipFile(file));
		}
		return zipFileInstances.get(file);
	}

	
	public void closeZipFiles() {
		
		// Try to close all zip files
		for (ZipFile zipFile : zipFileInstances.values()) {
			try {
				zipFile.close();
			} catch (Exception e) {
				logger.warn("could not close zip file");
			}
		}
		
		// Clear instance pool so that they are created again when needed next time
		zipFileInstances.clear();
		
	}
	
	@Override
	public void checkCompatibility(ContentLocation location) throws IllegalArgumentException {
		
		URL url = location.getUrl();

		// null url
		if (url == null) {
			throw new IllegalArgumentException("DataBean URL is null.");
		} 
		
		// protocol not "file"
		else if (!"file".equals(url.getProtocol())) {
			throw new IllegalArgumentException("Protocol of " + url.toString() + " is not \"file\".");
		} 
		
		// needs to be session file
		else if (!getZipFile(location).getName().endsWith("." + UserSession.SESSION_FILE_EXTENSION)) {
			throw new IllegalArgumentException("Not a session file.");
		}
		
		// needs to have non-empty reference
		else if (url.getRef() == null || url.getRef().length() == 0) {
			throw new IllegalArgumentException("Reference is null or empty.");
		}
	}

	
	private File getZipFile(ContentLocation location) {
		File zipFile;
		try {
			// remove fragment before converting to File
			URI beanURI = location.getUrl().toURI();
			URI zipURI = new URI(beanURI.getScheme(), beanURI.getSchemeSpecificPart(), null);
			zipFile = new File(zipURI);
		} catch (URISyntaxException use) {
			throw new IllegalArgumentException(location.getUrl() + " does not point to a file.");
		}
		return zipFile;
	}

	@Override
	public void markDeletable(ContentLocation location) {
		// entries are not deleted from zip files
	}

	@Override
	public OutputStream getOutputStream(ContentLocation location) throws IOException {
		throw new UnsupportedOperationException("zip content handler does not support output");
	}

	@Override
	public boolean isAccessible(ContentLocation location) {
		checkCompatibility(location);
		try {
			ZipFile zipFile = createZipFile(location);
			return zipFile.getEntry(location.getUrl().getRef()) != null;
		} catch (IOException e) {
			return false;
		}
	}

}
