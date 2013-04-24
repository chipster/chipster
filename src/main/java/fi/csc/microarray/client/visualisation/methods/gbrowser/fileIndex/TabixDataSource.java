package fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;

import org.broad.tribble.readers.TabixReader;

import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataSource;

/**
 * 
 * @author Petri Klemelä
 *
 */
public class TabixDataSource extends DataSource {

	private TabixReader reader;

    public TabixDataSource(URL tabixFile, URL tabixIndexFile) throws URISyntaxException, IOException {
    	//TODO use the provided index instead of guessing
        super(tabixFile, null);
        
        String fileString = null;
        
        if ("http".equals(tabixFile.getProtocol())) {
        	fileString = tabixFile.toExternalForm();
        } else {
        	fileString = (new File(tabixFile.toURI()).getPath()); //Translate '%20' to space character, required in Windows
        }
        
        this.reader = new TabixReader(fileString);

        // TODO check chromosome naming convention, see SAMDataSource
    }

	public TabixReader getReader() {
		return reader;
	}
}
