package fi.csc.microarray.filebroker;

import java.net.MalformedURLException;
import java.net.URL;
import java.util.Date;
import java.util.HashMap;
import java.util.UUID;

public class AuthorisedUrlRepository {
	
	private URL thisHost;

	private HashMap<URL, Date> repository = new HashMap<URL, Date>();  
	
	public AuthorisedUrlRepository(URL thisHost) {
		this.thisHost = thisHost;
		
	}

	public URL createAuthorisedUrl() throws MalformedURLException {
		
		// create unique url
		String filename = UUID.randomUUID().toString();
		URL newUrl;
		do {
			newUrl = new URL(thisHost.toString() + "/" + filename);
		} while (repository.containsKey(newUrl));
		
		// store it
		repository.put(newUrl, new Date());
		
		return newUrl;
	}
	

	public boolean isAuthorised(URL url) {
		for (URL key : repository.keySet()) {
			
		}
		
		throw new UnsupportedOperationException();
	}
	

}
