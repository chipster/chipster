package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URISyntaxException;
import java.net.URL;

public class DataUrl {

	private URL url;
	private String name;

	public DataUrl(URL data, String name) {
		this.url = data;
		this.name = name;
	}

	public String getName() {
		return name;
	}

	public InputStream getInputStream() throws IOException, URISyntaxException {

		//Assume local
		return new FileInputStream(new File(url.toURI()));
	}

	public File getLocalFile() throws IOException, URISyntaxException {
		//Assume local
		return new File(url.toURI());
	}

	public URL getUrl() throws IOException {
		return url;
	}
}