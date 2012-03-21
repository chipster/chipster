package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.FileNotFoundException;
import java.net.URISyntaxException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.TabixHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.TabixFile;

/**
 * Experimental data source for indexed Tabix files.
 * 
 * @author Petri Klemelä
 *
 */
public class TabixDataSource extends DataSource {
    
    private TabixFile tabixFile = null;
    
    
    /**
     * Generally we would like to have both data and index files,
     * because otherwise we could not access random locations.
     * 
     * @param file
     * @throws FileNotFoundException
     * @throws URISyntaxException 
     */
    public TabixDataSource(URL file) throws FileNotFoundException, URISyntaxException {
        super(file, TabixHandlerThread.class);
        tabixFile = new TabixFile(this.file);
    }
    
    public TabixFile getTabix() {
        return tabixFile;
    }
}