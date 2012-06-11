package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;

import org.broad.tribble.readers.TabixReader;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.GtfTabixHandlerThread;

/**
 * 
 * @author Petri Klemelä
 *
 */
public class GtfTabixDataSource extends DataSource {
	
	private TabixReader reader;

    /**
     * Generally we would like to have both data and index files,
     * because otherwise we could not access random locations.
     * 
     * @param file
     * @throws URISyntaxException 
     * @throws IOException 
     */
    public GtfTabixDataSource(URL tabixFile, URL tabixIndexFile) throws URISyntaxException, IOException {
    	//TODO use the provided index instead of guessing
        super(tabixFile, GtfTabixHandlerThread.class);
        
        String fileString = tabixFile.toString().replaceAll("file:/", "/");

        this.reader = new TabixReader(fileString);

        // TODO check chromosome naming convention, see SAMDataSource
    }

	public TabixReader getReader() {
		return reader;
	}
}
