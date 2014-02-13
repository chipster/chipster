package fi.csc.microarray.client.visualisation.methods.gbrowser.util;


public class UnsortedDataException extends GBrowserException {

	private String filename;

	public UnsortedDataException(String string, String filename) {
		super(string);
		this.filename = filename;
	}

	public String getFilename() {
		return filename;
	}
}
