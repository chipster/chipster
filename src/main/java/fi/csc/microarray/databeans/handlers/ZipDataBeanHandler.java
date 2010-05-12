package fi.csc.microarray.databeans.handlers;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.URI;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataBean.DataBeanType;

public class ZipDataBeanHandler extends DataBeanHandlerBase {

	public ZipDataBeanHandler() {
		super(DataBeanType.LOCAL_SESSION);
	}
	
	public long getContentLength(DataBean dataBean) throws IOException {
		checkCompatibility(dataBean);
		ZipFile zipFile = new ZipFile(getZipFile(dataBean));		
		ZipEntry zipEntry = zipFile.getEntry(dataBean.getContentUrl().getRef());
		return zipEntry.getSize();
	}

	public InputStream getInputStream(DataBean dataBean) throws IOException {
		checkCompatibility(dataBean);
		ZipFile zipFile = new ZipFile(getZipFile(dataBean));
		ZipEntry zipEntry = zipFile.getEntry(dataBean.getContentUrl().getRef());
		return zipFile.getInputStream(zipEntry);
	}

	protected void checkCompatibility(DataBean dataBean) throws IllegalArgumentException {
		super.checkCompatibility(dataBean);
		
		URL url = dataBean.getContentUrl();
		
		// null url
		if (url == null) {
			throw new IllegalArgumentException("DataBean URL is null.");
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
		else if (!getZipFile(dataBean).getName().endsWith(".cs")) {
			throw new IllegalArgumentException("Not a session file.");
		}
		
		// needs to have non-empty reference
		else if (url.getRef() == null || url.getRef().length() == 0) {
			throw new IllegalArgumentException("Reference is null or empty.");
		}
	}
	
	private File getZipFile(DataBean dataBean) {
		File zipFile;
		try {
			// remove fragment before converting to File
			URI beanURI = dataBean.getContentUrl().toURI();
			URI zipURI = new URI(beanURI.getScheme(), beanURI.getSchemeSpecificPart(), null);
			zipFile = new File(zipURI);
		} catch (URISyntaxException use) {
			throw new IllegalArgumentException(dataBean.getContentUrl() + " does not point to a file.");
		}
		return zipFile;
	}

	public void delete(DataBean dataBean) {
		// do nothing for now
	}

	public OutputStream getOutputStream(DataBean dataBean) throws IOException {
		// FIXME
		throw new UnsupportedOperationException("not supported yet");
	}
	
}
