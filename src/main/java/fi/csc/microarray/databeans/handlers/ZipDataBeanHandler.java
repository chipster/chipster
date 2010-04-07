package fi.csc.microarray.databeans.handlers;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import fi.csc.microarray.databeans.DataBean;

public class ZipDataBeanHandler implements DataBeanHandler {

	public long getContentLength(DataBean dataBean) throws IOException {
		checkURL(dataBean);
		ZipFile zipFile = new ZipFile(getFile(dataBean));		
		ZipEntry zipEntry = zipFile.getEntry(dataBean.getContentUrl().getRef());
		return zipEntry.getSize();
	}

	public InputStream getInputStream(DataBean dataBean) throws IOException {
		ZipFile zipFile = new ZipFile(getFile(dataBean));
		ZipEntry zipEntry = zipFile.getEntry(dataBean.getContentUrl().getRef());
		return zipFile.getInputStream(zipEntry);
	}

	private void checkURL(DataBean dataBean) throws IllegalArgumentException {
		URL url = dataBean.getContentUrl();
		
		// null url
		if (url == null) {
			throw new IllegalArgumentException("DataBean is null.");
		} 
		
		// protocol not "file"
		else if (!"file".equals(url.getProtocol())) {
			throw new IllegalArgumentException("Protocol of " + url.toString() + " is not \"file\".");
		} 
		
		// null or empty path
		else if (url.getPath() == null || url.getPath().length() == 0) {
			throw new IllegalArgumentException("Illegal path:" + url.toString());
		}
		
		// needs to be .cs
		else if (!getFile(dataBean).getName().endsWith(".cs")) {
			throw new IllegalArgumentException("Not a session file.");
		}
		
		// needs to have non-empty reference
		else if (url.getRef() == null || url.getRef().length() == 0) {
			throw new IllegalArgumentException("Reference is null or empty.");
		}
	}
	
	private File getFile(DataBean dataBean) {
		return new File(dataBean.getContentUrl().getPath());
	}

	
}
