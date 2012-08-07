package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.IndexedFastaHandlerThread;

/**
 * 
 * @author Petri Klemelä
 *
 */
public class IndexedFastaDataSource extends DataSource {

	private URL dataUrl;
    private URL index;

	public IndexedFastaDataSource(URL data, URL index) throws URISyntaxException, IOException {
        super(data, IndexedFastaHandlerThread.class);
        
        this.dataUrl = data;
        this.index = index;
    }
	
	public URL getIndex() {
		return index;
	}

	public URL getUrl() {
		return dataUrl;
	}
}
