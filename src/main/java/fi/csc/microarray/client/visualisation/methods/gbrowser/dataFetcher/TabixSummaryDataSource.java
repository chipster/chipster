package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.io.FileNotFoundException;
import java.net.URISyntaxException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;

/**
 * Data source for Tabix files.
 * 
 * @author naktinis, klemela
 *
 */
public class TabixSummaryDataSource extends DataSource {
    
    TabixFile tabixFile = null;
    
    
    public TabixSummaryDataSource(URL url) throws FileNotFoundException, URISyntaxException, GBrowserException {
        // TODO Support URLs            
        super(url, TabixSummaryHandlerThread.class);
        throw new GBrowserException("Reading Tabix files from URL not yet supported");
    }

    public TabixFile getTabix() {
        return tabixFile;
    }
}
