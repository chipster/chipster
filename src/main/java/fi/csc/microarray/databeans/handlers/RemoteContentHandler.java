package fi.csc.microarray.databeans.handlers;

import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.Calendar;
import java.util.Date;
import java.util.LinkedList;

import fi.csc.microarray.databeans.DataManager.ContentLocation;
import fi.csc.microarray.util.KeyAndTrustManager;
import fi.csc.microarray.util.UrlTransferUtil;

public class RemoteContentHandler implements ContentHandler {
	
	private static int BLACKLIST_LIFETIME_MILLISECONDS = 30*UrlTransferUtil.HTTP_TIMEOUT_MILLISECONDS;

	private LinkedList<BlacklistEntry> hostBlacklist = new LinkedList<BlacklistEntry>();
	
	private static class BlacklistEntry {

		public String host;
		public int port;
		public Date created = new Date();
		
		public BlacklistEntry(String host, int port) {
			this.host = host;
			this.port = port;
		}
		
		@Override
		public boolean equals(Object other) {
			
			if (other instanceof BlacklistEntry) {
				BlacklistEntry otherEntry = (BlacklistEntry)other;
				
				return this.host.equals(otherEntry.host) && this.port == otherEntry.port;
			} else {
				return false;
			}
		}
		
	}
	
	@Override
	public InputStream getInputStream(ContentLocation location) throws IOException {
		checkCompatibility(location);
		HttpURLConnection connection = (HttpURLConnection)location.getUrl().openConnection();
		// filebroker isn't a ContentLocation, so use CA certs
		KeyAndTrustManager.configureForCACertificates(connection);
		return connection.getInputStream();
	}

	@Override
	public OutputStream getOutputStream(ContentLocation location) throws IOException {
		throw new UnsupportedOperationException("remote content handler does not support output");
	}
	
	@Override
	public Long getContentLength(ContentLocation location) throws IOException {
		checkCompatibility(location);
		// filebroker isn't a ContentLocation, so use CA certs
		return UrlTransferUtil.getContentLength(location.getUrl(), false);
	}

	/**
	 * Deletes files that have been created for this bean. Only deletes temporary local files,
	 * not user files etc.
	 */
	@Override
	public void markDeletable(ContentLocation location) {
		// remote data is not currently deleted
	}
	
	@Override
	public void checkCompatibility(ContentLocation location) throws IllegalArgumentException {
		
		URL url = location.getUrl();
		
		// null url
		if (url == null) {
			throw new IllegalArgumentException("url is null.");
		} 

		// protocol not HTTP compatible
		else if (!("http".equals(url.getProtocol()) || "https".equals(url.getProtocol()))) {
			throw new IllegalArgumentException("Protocol of " + url.toString() + " is not HTTP compatible.");
		}
		
		// null or empty path
		else if (url.getPath() == null || url.getPath().length() == 0) {
			throw new IllegalArgumentException("Illegal path: " + url.toString());
		} 


	}

	@Override
	public boolean isAccessible(ContentLocation location) {
		checkCompatibility(location);
		URL url = location.getUrl();

		try {
			
			if (isBlacklisted(url.getHost(), url.getPort())) {
				return false; // do not check blacklisted hosts
				
			} else {
				// filebroker isn't a ContentLocation, so use CA certs
				return UrlTransferUtil.isAccessible(url, false);
			}
			
		} catch (IOException e) {
			// timeout or other IO problem when connecting to server => don't try again
			blacklist(url.getHost(), url.getPort());
			return false;
		}
	}
	
	private boolean isBlacklisted(String host, int port) {
		
		// first prune outdated
		Calendar c = Calendar.getInstance();
		c.add(Calendar.MILLISECOND, -BLACKLIST_LIFETIME_MILLISECONDS);
		Date lastValid = c.getTime();
		LinkedList<BlacklistEntry> outdated = new LinkedList<RemoteContentHandler.BlacklistEntry>();
		for (BlacklistEntry entry : hostBlacklist) {
			if (entry.created.before(lastValid)) {
				outdated.add(entry);
			}
		}
		hostBlacklist.removeAll(outdated);
		
		// check if listed
		return hostBlacklist.contains(new BlacklistEntry(host, port));
	}
	

	private void blacklist(String host, int port) {
		hostBlacklist.add(new BlacklistEntry(host, port));
	}
}
