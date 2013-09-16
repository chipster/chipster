package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.MalformedURLException;
import java.net.URISyntaxException;
import java.net.URL;

/**
 * Genome browser shouldn't care whether the data is a local file or an url address. This class
 * is a replacement for File and URL classes to enable transparent use of data regardless of its
 * type. 
 * 
 * @author klemela
 */
public class DataUrl {

	private URL url;
	private String name;

	public DataUrl(URL data, String name) {
		this.url = data;
		this.name = name;
	}

	public DataUrl(File file) throws MalformedURLException {
		this(file.toURI().toURL(), file.getName());
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

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result + ((name == null) ? 0 : name.hashCode());
		result = prime * result + ((url == null) ? 0 : url.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (obj == null) {
			return false;
		}
		if (!(obj instanceof DataUrl)) {
			return false;
		}
		DataUrl other = (DataUrl) obj;
		if (name == null) {
			if (other.name != null) {
				return false;
			}
		} else if (!name.equals(other.name)) {
			return false;
		}
		if (url == null) {
			if (other.url != null) {
				return false;
			}
		} else if (!url.equals(other.url)) {
			return false;
		}
		return true;
	}		
}