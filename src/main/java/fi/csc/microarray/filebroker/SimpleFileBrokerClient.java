package fi.csc.microarray.filebroker;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;
import java.util.LinkedList;

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

	@Override
	public URL addFile(FileBrokerArea area, InputStream content, long contentLength, CopyProgressListener progressListener) throws FileBrokerException, JMSException, IOException {
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
	public URL getPublicUrl() throws MalformedURLException {
		return new URL(DirectoryLayout.getInstance().getConfiguration().getString("messaging", "public-files-url"));
	}

	@Override
	public void getFile(File file, URL inputUrl) throws IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public URL addFile(FileBrokerArea area, File file, CopyProgressListener progressListener) throws FileBrokerException, JMSException, IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean requestDiskSpace(long size) {
		throw new UnsupportedOperationException();
	}

	@Override
	public URL moveFileToStorage(URL url, long contentLength) {
		throw new UnsupportedOperationException();
	}

	@Override
	public String[][] listRemoteSessions() throws JMSException {
		throw new UnsupportedOperationException();
	}

	@Override
	public void removeRemoteSession(URL sessionURL) throws JMSException {
		throw new UnsupportedOperationException();
	}

	@Override
	public void saveRemoteSession(String name, URL sessionURL, LinkedList<URL> dataUrls)
			throws JMSException {
		throw new UnsupportedOperationException();
	}

	@Override
	public URL addSessionFile() throws JMSException, FileBrokerException {
		throw new UnsupportedOperationException();
	}

}
