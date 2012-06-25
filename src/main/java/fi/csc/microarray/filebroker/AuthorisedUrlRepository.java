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

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.security.CryptoKey;

/**
 * Stores URL's that have been authorised to be used for PUT'ing content into. Objects
 * of this class are thread safe.
 * 
 * @author Aleksi Kallio
 *
 */
public class AuthorisedUrlRepository {

	/**
	 * The time after which URL is removed from the repository, 
	 * i.e., PUT requests are to be rejected to that URL.
	 */
	private static final int URL_LIFETIME_MINUTES = 10;
	private static final String COMPRESSION_SUFFIX = ".compressed";

	private HashMap<URL, Date> repository = new HashMap<URL, Date>();  
	private Lock repositoryLock = new ReentrantLock();

	private String host;
	private int port;
	private String cachePath;
	
	public AuthorisedUrlRepository(String host, int port) {
		this.host = host;
		this.port = port;
		
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
		cachePath = configuration.getString("filebroker", "cache-path");
	}

	/**
	 * Creates new URL and adds it to repository, where it has a 
	 * limited lifetime.
	 * 
	 *  @see #URL_LIFETIME_MINUTES
	 */
	public URL createAuthorisedUrl(boolean useCompression) throws MalformedURLException {

		URL newUrl;

		String compressionSuffix = "";
		if (useCompression) {
			compressionSuffix = COMPRESSION_SUFFIX;
		}

		repositoryLock.lock();
		try {
			// create url that does not exist in the repository
			do {
				String filename = CryptoKey.generateRandom();
				newUrl = new URL(host + ":" + port + "/" + cachePath + "/" + filename + compressionSuffix);
				
			} while (repository.containsKey(newUrl));

			// store it
			repository.put(newUrl, new Date());
			
		} finally {
			repositoryLock.unlock();
		}
		
		return newUrl;
	}
	
	/**
	 * Checks if repository contains valid (not outdated) copy of the URL. 
	 */
	public boolean isAuthorised(URL url) {
	
		boolean contains = false;

		repositoryLock.lock();
		try {
			// prune non-valid dates first
			Iterator<URL> keyIterator = repository.keySet().iterator(); // use iterator because we are removing
			while (keyIterator.hasNext()) {
				URL key = keyIterator.next();
				if (!isDateValid(repository.get(key))) {
					keyIterator.remove();
				}
			}

			// check if url exists in up-to-date repository
			contains = repository.containsKey(url);

		} finally {
			repositoryLock.unlock();
		}

		return contains;
	}
	
	public boolean checkFilenameSyntax(String filename) {
		String fileNameToCheck = filename;
		if (filename.endsWith(COMPRESSION_SUFFIX)) {
			fileNameToCheck = filename.substring(0, filename.length() - COMPRESSION_SUFFIX.length());
		}
		
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


