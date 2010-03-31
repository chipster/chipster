package fi.csc.microarray.databeans.handlers;

import java.io.BufferedInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;

import fi.csc.microarray.databeans.DataBean;

public class LocalFileDataBeanHandler implements DataBeanHandler {

	public InputStream getInputStream(DataBean dataBean) throws FileNotFoundException {
		checkURL(dataBean);
		return new BufferedInputStream(new FileInputStream(getFile(dataBean)));
	}

	public long getContentLength(DataBean dataBean) {
		checkURL(dataBean);
		return getFile(dataBean).length();
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
	}
	
	private File getFile(DataBean dataBean) {
		return new File(dataBean.getContentUrl().getPath());
	}

	public static void main(String[] args) throws MalformedURLException {
		URL url1 = new URL("file:///home/hupponen/neppi.txt");
		
		URL url = new URL(url1, "#jepjep");
		System.out.println(url.getProtocol());
		System.out.println(url.getPath());
		System.out.println(url.getFile());
		System.out.println(url.getUserInfo());
		System.out.println(url.getRef());
	}

}
