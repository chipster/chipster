package fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex;

import java.io.File;
import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;

import org.broad.tribble.readers.TabixReader;
import org.broad.tribble.readers.TabixReader.Iterator;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.SamBamUtils;

/**
 * 
 * @author Petri Klemelä
 *
 */
public class TabixDataSource extends DataSource {

	private TabixReader reader;

    public TabixDataSource(DataUrl repeat, DataUrl repeatIndex) throws URISyntaxException, IOException {
    	//TODO use the provided index instead of guessing
        super(repeat);
        
        URL repeatUrl = repeat.getUrl();
        
        String fileString = null;
        
        if ("http".equals(repeatUrl.getProtocol())) {
        	fileString = repeatUrl.toExternalForm();
        } else {
        	fileString = (new File(repeatUrl.toURI()).getPath()); //Translate '%20' to space character, required in Windows
        }
        
        this.reader = new TabixReader(fileString);

        // TODO initialize chromosome name unnormaliser (see for example BamDataSource), 
        //However, it isn't possible to get list of chromosomes from TabixReader. We could read that list 
        //from the indexFile. This is not strictly necessary at the moment, because all our gtf 
        //and repeat files obey Ensembl naming convention.
    }
    
	public void clean() {
		SamBamUtils.closeIfPossible(reader);
	}
	
	public Iterator getTabixIterator(Region request) {
		String chromosome = request.start.chr.toNormalisedString();

		//limit to integer range
		int start = (int) Math.min(Integer.MAX_VALUE, request.start.bp);
		int end = (int) Math.min(Integer.MAX_VALUE, request.end.bp);
		
		//Extend area to be able to draw regions that start before the left screen edge, but don't go over MAX_VALUE, or below 1
		//TODO Be more clever to avoid getting so much useless data
		int EXTRA = 500000; //O,5M should be enought for the longest human introns http://www.bioinfo.de/isb/2004040032/
		
		start = (int) Math.max((long)start - EXTRA, 1);
		end = (int) Math.min((long)end + EXTRA, Integer.MAX_VALUE);

		//Check that region is below max bin size of Tabix
		int MAX_BIN_SIZE = 512*1024*1024 - 2;

		start = (int) Math.min(MAX_BIN_SIZE, start);
		end = (int) Math.min(MAX_BIN_SIZE, end);

		start = (int) Math.max(1, start);
		end = (int) Math.max(1, end);

		String queryRegion = chromosome + ":" + start + "-" + end;

		TabixReader.Iterator iter = null;
		
		try {
			iter = reader.query(queryRegion);
		} catch (ArrayIndexOutOfBoundsException e) {
			//No such chromosome
		}
			
		return iter;
	}
}
