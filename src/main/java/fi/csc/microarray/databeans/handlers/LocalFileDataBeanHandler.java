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
import java.net.URL;

import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.databeans.DataBean.StorageMethod;

public class LocalFileDataBeanHandler extends DataBeanHandlerBase {

	public LocalFileDataBeanHandler(DataManager dataManager) {
		super(dataManager, StorageMethod.LOCAL_USER, StorageMethod.LOCAL_TEMP);
	}
	
	
	public InputStream getInputStream(DataBean dataBean) throws FileNotFoundException {
		checkCompatibility(dataBean);
		return new BufferedInputStream(new FileInputStream(getFile(dataBean)));
	}

	public OutputStream getOutputStream(DataBean dataBean) throws IOException {
		checkCompatibility(dataBean);
		return new BufferedOutputStream(new FileOutputStream(getFile(dataBean)));
	}
	
	public long getContentLength(DataBean dataBean) {
		checkCompatibility(dataBean);
		return getFile(dataBean).length();
	}

	/**
	 * Only delete temporary files, never user files.
	 */
	public void delete(DataBean dataBean) {
		if (dataBean.getStorageMethod().equals(DataBean.StorageMethod.LOCAL_TEMP)) {

			checkCompatibility(dataBean);
			File file = getFile(dataBean);
			file.delete();
		}
	}
	
	
	protected void checkCompatibility(DataBean dataBean) throws IllegalArgumentException {
		super.checkCompatibility(dataBean);

		URL url = dataBean.getContentUrl();
		
		// null url
		if (url == null) {
			throw new IllegalArgumentException("DataBean url is null.");
		} 
		
		// protocol not "file"
		else if (!"file".equals(url.getProtocol())) {
			throw new IllegalArgumentException("Protocol of " + url.toString() + " is not \"file\".");
		} 
		
		// null or empty path
		else if (url.getPath() == null || url.getPath().length() == 0) {
			throw new IllegalArgumentException("Illegal path:" + url.toString());
		} 
	}
	
	public File getFile(DataBean dataBean) {
		File file;
		try {
			file = new File(dataBean.getContentUrl().toURI());
		} catch (URISyntaxException use) {
			throw new IllegalArgumentException(dataBean.getContentUrl() + " does not point to a file.");
		}
		return file;
	}

}
