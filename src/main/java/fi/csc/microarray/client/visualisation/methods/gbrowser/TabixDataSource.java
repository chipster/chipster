package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.TabixFile;
import fi.csc.microarray.exception.MicroarrayException;

/**
 * Data source for Tabix files.
 * 
 * @author naktinis, klemela
 *
 */
public class TabixDataSource extends DataSource {
    
    TabixFile tabixFile = null;
    
    
    public TabixDataSource(URL url) throws FileNotFoundException, MicroarrayException {
        // TODO Support URLs    	
        super(url);
        throw new MicroarrayException("Reading Tabix files from URL not yet supported");
    }

    /**
     * Generally we would like to have both data and index files,
     * because otherwise we could not access random locations.
     * 
     * @param file
     * @throws FileNotFoundException
     */
    public TabixDataSource(File file) throws FileNotFoundException {
        super(file);
        tabixFile = new TabixFile(file);
    }
    
    public TabixFile getTabix() {
        return tabixFile;
    }
}