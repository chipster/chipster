package fi.csc.microarray.filebroker;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.LinkedList;
import java.util.List;

import javax.jms.JMSException;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.util.IOUtils.CopyProgressListener;

/**
 * Simple file broker client for the standalone mode.
 * 
 * Only supports getPublicUrl() for now.
 * 
 * @author hupponen
 *
 */
public class SimpleFileBrokerClient implements FileBrokerClient {
	
	private static final String PUBLIC_FILES = "public-files.txt";

	@Override
	public URL addFile(InputStream content, long contentLength, CopyProgressListener progressListener) throws FileBrokerException, JMSException, IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean checkFile(URL url, long contentLength) {
		throw new UnsupportedOperationException();	}

	@Override
	public InputStream getFile(URL url) throws IOException {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public List<URL> getPublicFiles() throws JMSException, MalformedURLException {
		
		String publicRoot = DirectoryLayout.getInstance().getConfiguration().getString("messaging", "public-files-url") + "/";
		URL filesListing = new URL(publicRoot + PUBLIC_FILES);

		List<URL> list = new LinkedList<URL>();
		
		try {
			
			BufferedReader reader = new BufferedReader(new InputStreamReader(filesListing.openStream()));

			String line;
			while ((line = reader.readLine()) != null) {
				list.add(new URL(publicRoot + line));
			}
			
		} catch (IOException e) {
			throw new IllegalStateException("Unable to read public file list from server", e);
		}
			
		return list;		
	}

	@Override
	public URL getPublicUrl() throws MalformedURLException {
		return new URL(DirectoryLayout.getInstance().getConfiguration().getString("messaging", "public-files-url"));
	}

	@Override
	public void getFile(File file, URL inputUrl) throws IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public URL addFile(File file, CopyProgressListener progressListener) throws FileBrokerException, JMSException, IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean requestDiskSpace(long size) {
		throw new UnsupportedOperationException();
	}

}
