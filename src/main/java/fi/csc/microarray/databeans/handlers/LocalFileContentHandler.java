package fi.csc.microarray.databeans.handlers;

import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URISyntaxException;

import fi.csc.microarray.databeans.DataBean.ContentLocation;
import fi.csc.microarray.databeans.DataManager.StorageMethod;

public class LocalFileContentHandler implements ContentHandler {
	
	public InputStream getInputStream(ContentLocation location) throws FileNotFoundException {
		return new BufferedInputStream(new FileInputStream(getFile(location)));
	}

	public OutputStream getOutputStream(ContentLocation location) throws IOException {
		return new BufferedOutputStream(new FileOutputStream(getFile(location)));
	}
	
	public long getContentLength(ContentLocation location) {
		return getFile(location).length();
	}

	/**
	 * Deletes files that have been created for this bean. Only deletes temporary local files,
	 * not user files etc.
	 */
	public void canBeDeleted(ContentLocation location) {
		
		if (location.getMethod() == StorageMethod.LOCAL_TEMP) {
			File file = getFile(location);
			file.delete();
		}
	}
	
	
//	protected void checkCompatibility(DataBean dataBean) throws IllegalArgumentException {
//
//		URL url = dataBean.getLocalUrl();
//		
//		// null url
//		if (url == null) {
//			throw new IllegalArgumentException("DataBean url is null.");
//		} 
//		
//		// protocol not "file"
//		else if (!"file".equals(url.getProtocol())) {
//			throw new IllegalArgumentException("Protocol of " + url.toString() + " is not \"file\".");
//		} 
//		
//		// null or empty path
//		else if (url.getPath() == null || url.getPath().length() == 0) {
//			throw new IllegalArgumentException("Illegal path:" + url.toString());
//		} 
//	}
	
	public File getFile(ContentLocation location) {
		File file;
		try {
			file = new File(location.getUrl().toURI());
		} catch (URISyntaxException use) {
			throw new IllegalArgumentException(location.getUrl() + " does not point to a file.");
		}
		return file;
	}

}
