package fi.csc.microarray.client.visualisation.methods.gbrowser.util;

public class GBrowserException extends Exception {

	public GBrowserException(String string) {
		super(string);
	}

	public GBrowserException(String string, Exception e) {
		super(string, e);
	}
}
