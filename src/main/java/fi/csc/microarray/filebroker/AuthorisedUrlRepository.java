package fi.csc.microarray.filebroker;

import java.net.MalformedURLException;
import java.net.URL;
import java.util.Calendar;
import java.util.Date;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.UUID;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

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

	private URL thisHost;
	private HashMap<URL, Date> repository = new HashMap<URL, Date>();  
	private Lock repositoryLock = new ReentrantLock();
	
	public AuthorisedUrlRepository(URL thisHost) {
		this.thisHost = thisHost;
		
	}

	/**
	 * Creates new URL and adds it to repository, where it has a 
	 * limited lifetime.
	 * 
	 *  @see #URL_LIFETIME_MINUTES
	 */
	public URL createAuthorisedUrl() throws MalformedURLException {

		URL newUrl;

		repositoryLock.lock();
		try {
			// create url that does not exist in the repository
			do {
				String filename = UUID.randomUUID().toString();
				newUrl = new URL(thisHost.toString() + "/" + filename);
				
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
			for (URL key : repository.keySet()) {
				if (!isDateValid(repository.get(key))) {
					repository.remove(key);
				}
			}

			// check if url exists in up-to-date repository
			contains = repository.containsKey(url);

		} finally {
			repositoryLock.unlock();
		}

		return contains;
	}
	
	private boolean isDateValid(Date date) {
		Calendar cal = new GregorianCalendar();
		// go back from current time 
		cal.add(Calendar.MINUTE, -URL_LIFETIME_MINUTES); 		
		return cal.getTime().before(date);
	}
}
