package fi.csc.microarray.filebroker;

import java.net.MalformedURLException;
import java.net.URL;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.Iterator;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import org.apache.log4j.Logger;

import fi.csc.microarray.filebroker.FileBrokerClient.FileBrokerArea;
import fi.csc.microarray.security.CryptoKey;

/**
 * Stores URL's that have been authorised to be used for PUT'ing content into. Objects
 * of this class are thread safe.
 * 
 * @author Aleksi Kallio
 *
 */
public class AuthorisedUrlRepository {
	
	public static class Authorisation {
		private Date created;
		private long fileSize;
		
		public Authorisation(Date created, long fileSize) {
			this.created = created;
			this.fileSize = fileSize;
		}
		
		public Authorisation(long fileSize) {
			this(new Date(), fileSize);
		}
		
		public Authorisation(Authorisation o) {
			this(o.created, o.fileSize);
		}

		public Date getCreated() {
			return created;
		}
		public void setCreated(Date created) {
			this.created = created;
		}

		public long getFileSize() {
			return fileSize;
		}

		public void setFileSize(long fileSize) {
			this.fileSize = fileSize;
		}
	}
	
	private static Logger logger = Logger.getLogger(AuthorisedUrlRepository.class);

	/**
	 * The time after which URL is removed from the repository, 
	 * i.e., PUT requests are to be rejected to that URL.
	 */
	private static final int URL_LIFETIME_MINUTES = 10;
	private static final String COMPRESSION_SUFFIX = ".compressed";

	private HashMap<URL, Authorisation> repository = new HashMap<>();  
	private Lock repositoryLock = new ReentrantLock();

	private String host;
	private int port;
	private String cachePath;
	private String storagePath;
	
	public AuthorisedUrlRepository(String host, int port, String cachePath, String storagePath) {
		this.host = host;
		this.port = port;
		this.cachePath = cachePath;
		this.storagePath = storagePath;
	}

	/**
	 * Creates new URL and adds it to repository, where it has a 
	 * limited lifetime.
	 * @param area 
	 * @param bytes 
	 * 
	 *  @see #URL_LIFETIME_MINUTES
	 */
	public URL createAuthorisedUrl(String fileId, boolean useCompression, FileBrokerArea area, long bytes) throws Exception {

		URL newUrl;

		String compressionSuffix = "";
		if (useCompression) {
			compressionSuffix = COMPRESSION_SUFFIX;
		}

		repositoryLock.lock();
		try {
			// create url
			if (area == FileBrokerArea.STORAGE) {
				newUrl = constructStorageURL(fileId, compressionSuffix);
			} else {
				newUrl = constructCacheURL(fileId, compressionSuffix);
			}
				
			if (repository.containsKey(newUrl)) {
				//two clients tried to upload the same file at the same time
				//allow and handle collision in RestServlet
				logger.debug("repository contains the url already, timestamp is renewed");
			}

			// store it
			repository.put(newUrl, new Authorisation(bytes));
			
		} finally {
			repositoryLock.unlock();
		}
		
		return newUrl;
	}
	
	public URL constructStorageURL(String filename, String compressionSuffix) throws MalformedURLException {
		return new URL(host + ":" + port + "/" + storagePath + "/" + filename + compressionSuffix);
	}

	public URL constructCacheURL(String filename, String compressionSuffix) throws MalformedURLException {
		return new URL(host + ":" + port + "/" + cachePath + "/" + filename + compressionSuffix);
	}

	/**
	 * Checks if repository contains valid (not outdated) copy of the URL.
	 * @param url
	 * @return Authorisation object the URL is valid, otherwise null
	 */
	public Authorisation getAuthorisation(URL url) {
	
		Authorisation authorisation = null;

		repositoryLock.lock();
		try {
			// prune non-valid dates first
			Iterator<URL> keyIterator = repository.keySet().iterator(); // use iterator because we are removing
			while (keyIterator.hasNext()) {
				URL key = keyIterator.next();
				if (!isDateValid(repository.get(key).getCreated())) {
					keyIterator.remove();
				}
			}

			// check if url exists in up-to-date repository
			authorisation = repository.get(url);

		} finally {
			repositoryLock.unlock();
		}

		if (authorisation != null) {
			//clone
			return new Authorisation(authorisation);
		} else {
			return null;
		}
	}
	
	
	public static String stripCompressionSuffix(String filename) {
		if (filename.endsWith(COMPRESSION_SUFFIX)) {
			return filename.substring(0, filename.length() - COMPRESSION_SUFFIX.length());
		} else {
			return filename;
		}

	}
	public static boolean checkFilenameSyntax(String filename) {
		String fileNameToCheck = stripCompressionSuffix(filename);
		return CryptoKey.validateKeySyntax(fileNameToCheck);
	}
	
	private boolean isDateValid(Date date) {
		Calendar cal = new GregorianCalendar();
		// go back from current time 
		cal.add(Calendar.MINUTE, -URL_LIFETIME_MINUTES); 		
		return cal.getTime().before(date);
	}

	public String getRootUrl() {
		return host + ":" + port;
	}
	
}


