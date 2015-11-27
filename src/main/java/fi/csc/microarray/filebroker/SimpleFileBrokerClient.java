package fi.csc.microarray.filebroker;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.MalformedURLException;
import java.net.URL;
import java.net.URLConnection;
import java.util.LinkedList;
import java.util.List;
import java.util.UUID;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.messaging.admin.StorageAdminAPI.StorageEntryMessageListener;
import fi.csc.microarray.util.IOUtils.CopyProgressListener;
import fi.csc.microarray.util.KeyAndTrustManager;

/**
 * Simple file broker client for the standalone mode.
 * 
 * Only supports getPublicFiles() for now.
 * 
 * @author hupponen
 *
 */
public class SimpleFileBrokerClient implements FileBrokerClient {

	private static final String PUBLIC_FILES = "public-files.txt";

	@Override
	public String addFile(String dataId, FileBrokerArea area, InputStream content, long contentLength, CopyProgressListener progressListener) throws FileBrokerException, IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public ChecksumInputStream getInputStream(String dataId) throws IOException {
		throw new UnsupportedOperationException();
	}
	
	@Override
	public List<URL> getPublicFiles() throws FileBrokerException, MalformedURLException {
		
		String publicRoot = DirectoryLayout.getInstance().getConfiguration().getString("messaging", "public-files-url") + "/";
		URL filesListing = new URL(publicRoot + PUBLIC_FILES);

		List<URL> list = new LinkedList<URL>();
		
		try {
			URLConnection connection = filesListing.openConnection();
			KeyAndTrustManager.configureForChipsterCertificate(connection);
			BufferedReader reader = new BufferedReader(new InputStreamReader(connection.getInputStream()));

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
	public void getFile(UUID sessionId, String dataId, File file) throws IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public void addFile(UUID sessionId, String dataId, FileBrokerArea area, File file, CopyProgressListener progressListener) throws FileBrokerException, FileBrokerException, IOException {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean requestDiskSpace(long size) {
		throw new UnsupportedOperationException();
	}

	@Override
	public List<DbSession> listRemoteSessions() throws FileBrokerException {
		throw new UnsupportedOperationException();
	}

	@Override
	public void removeRemoteSession(String dataId) throws FileBrokerException {
		throw new UnsupportedOperationException();
	}

	@Override
	public void saveRemoteSession(String name, String dataId, LinkedList<String> dataIds)
			throws FileBrokerException {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean isAvailable(String dataId, Long size, String checksum, FileBrokerArea area) {
		throw new UnsupportedOperationException();
	}

	@Override
	public boolean moveFromCacheToStorage(String dataId) {
		throw new UnsupportedOperationException();
	}

	@Override
	public List<DbSession> listPublicRemoteSessions() throws FileBrokerException {
		throw new UnsupportedOperationException();
	}

	@Override
	public String getExternalURL(String dataId) throws FileBrokerException,
			FileBrokerException {
		throw new UnsupportedOperationException();
	}

	@Override
	public Long getContentLength(String dataId) throws IOException,
	FileBrokerException, FileBrokerException {
		throw new UnsupportedOperationException();
	}

	@Override
	public StorageEntryMessageListener getStorageUsage() throws FileBrokerException,
			InterruptedException {
		throw new UnsupportedOperationException();
	}
}
