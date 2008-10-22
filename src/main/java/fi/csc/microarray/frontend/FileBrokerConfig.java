package fi.csc.microarray.frontend;

import java.net.MalformedURLException;
import java.net.URL;
import java.util.Random;

import fi.csc.microarray.MicroarrayConfiguration;

public class FileBrokerConfig {
	
	private static final String URL_SEPARATOR = "/";
	
	private static Random randGenerator = new Random();
	private static String[] FILESERVER_URLS;

	static {
		FILESERVER_URLS = MicroarrayConfiguration.getValues("messaging", "filebroker_urls");
		
		// add dir separators if they are missing
		for (int i = 0; i < FILESERVER_URLS.length; i++) {
			if (!FILESERVER_URLS[i].endsWith(URL_SEPARATOR)) {
				FILESERVER_URLS[i] = FILESERVER_URLS[i] + URL_SEPARATOR;
			}
		}
	}
	
	/**
	 * Selects one configured file broker by random and generates URL for given file
	 * to that broker.
	 */
	public static URL generateUrlToSomeFileBroker(String filename) throws MalformedURLException {
		int i = randGenerator.nextInt(FILESERVER_URLS.length);
		return new URL(FILESERVER_URLS[i] + filename);
	}

	/**
	 * Returns port for this file broker instance. Should be used only by filebrokers.
	 */
	public static int getPort() {
		return Integer.parseInt(MicroarrayConfiguration.getValue("filebroker", "port"));
	}
}
