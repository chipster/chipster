package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.net.URL;

public class SingleFile extends FileConfiguration {

	private URL file;

	public SingleFile (URL file) {
		this.file = file;
	}
	
	public URL getFile() {
		return file;
	}
}
