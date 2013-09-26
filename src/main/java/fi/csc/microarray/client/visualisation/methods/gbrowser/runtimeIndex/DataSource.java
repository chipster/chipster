package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;

/**
 * <p>One source of genomic content. Abstraction hides the physical data source, that can be either 
 * a file or a URL accessed via HTTP.</p>
 *
 * <p>TODO: add comparison methods so that when DataSource is used as a key
 * in a hash map, data sources with identical urls or file paths would be
 * considered the same.</p>
 * 
 * @author Petri Klemelä
 *
 */
public abstract class DataSource {

	protected DataUrl dataUrl;
	protected File file = null;
	protected URL url = null;
	protected String name;

	public DataSource(DataUrl dataUrl) throws URISyntaxException, IOException {
		
		this.dataUrl = dataUrl;
		this.url = dataUrl.getUrl();
		
		if (url != null) {
			if ("file".equals(url.getProtocol())) {
				file = new File(url.toURI());
			}
			this.name = url.toString(); 
		}
	}
	
	@Override
	public String toString() {
		return name;
	}

	public DataUrl getDataUrl() {
		return dataUrl;
	}
}
