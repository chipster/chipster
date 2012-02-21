package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.File;
import java.io.FileNotFoundException;

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
     */
    public TabixDataSource(File file) throws FileNotFoundException {
        super(file);
        tabixFile = new TabixFile(file);
    }
    
    public TabixFile getTabix() {
        return tabixFile;
    }
}