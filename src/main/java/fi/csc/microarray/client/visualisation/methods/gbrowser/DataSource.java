package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.MalformedURLException;
import java.net.URL;

/**
 * One source of genomic content, typically corresponding to one track. Abstraction hides
 * the physical data source, that can be either a file or a URL accessed via HTTP.
 *
 */
public abstract class DataSource {

	protected File file = null;
	protected URL url = null;
	protected String name;

	public DataSource(URL url) throws FileNotFoundException {
		this.url = url;
		this.name = url.toString(); 
	}

	public DataSource(File file) throws FileNotFoundException {
		this.file = file;
		this.name = file.toString();
	}
	
	public DataSource(URL urlRoot, String path)
	        throws FileNotFoundException, MalformedURLException {
		this(new URL(urlRoot.toString() + "/" + path));
	}

	public DataSource(File fileRoot, String path)
	        throws FileNotFoundException, MalformedURLException {
		this(new File(fileRoot, path));
	}
	
	@Override
	public String toString() {
		return name;
	}
}
