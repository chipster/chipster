package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.File;
import java.io.FileNotFoundException;
import java.net.URL;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.SAMFile;

/**
 * Data source for SAM files.
 * 
 * @author naktinis
 *
 */
public class SAMDataSource extends DataSource {
    
    SAMFile samFile = null;
    
    
    public SAMDataSource(URL url) throws FileNotFoundException {
        // TODO Support URLs
        super(url);
    }

    /**
     * Generally we would like to have both data and index files,
     * because otherwise we could not access random locations.
     * 
     * @param file
     * @throws FileNotFoundException
     */
    public SAMDataSource(File file, File indexFile) throws FileNotFoundException {
        super(file);
        samFile = new SAMFile(file, indexFile);
    }
    
    public SAMFile getSAM() {
        return samFile;
    }
}
