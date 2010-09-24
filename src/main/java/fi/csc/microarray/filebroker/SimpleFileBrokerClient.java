package fi.csc.microarray.filebroker;

import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URL;

import javax.jms.JMSException;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.util.IOUtils.CopyProgressListener;

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

}
