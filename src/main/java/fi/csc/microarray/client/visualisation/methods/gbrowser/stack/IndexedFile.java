package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.net.URL;

public class IndexedFile extends FileConfiguration{

	private URL data;
	private URL index;

	public IndexedFile(URL data, URL index) {
		this.data = data;
		this.index = index;
	}
	
	public URL getData() {
		return data;
	}
	
	public URL getIndex() {
		return index;
	}
}
