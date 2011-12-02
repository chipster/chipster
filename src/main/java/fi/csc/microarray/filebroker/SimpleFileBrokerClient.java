package fi.csc.microarray.filebroker;

import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;

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
	public URL addFile(InputStream content, CopyProgressListener progressListener) throws FileBrokerException, JMSException, IOException {
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
	public void getLocalFile(File file, URL inputUrl) throws IOException {
		throw new UnsupportedOperationException();
	}

}
