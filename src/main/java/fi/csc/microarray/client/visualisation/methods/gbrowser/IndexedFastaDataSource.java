package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.IOException;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Iterator;

import net.sf.picard.reference.ChipsterIndexedFastaSequenceFile;
import net.sf.samtools.SAMSequenceRecord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.IndexedFastaHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ChromosomeNameUnnormaliser;

/**
 * 
 * @author Petri Klemelä
 *
 */
public class IndexedFastaDataSource extends DataSource {
	
	private ChromosomeNameUnnormaliser chromosomeNameUnnormaliser = ChromosomeNameUnnormaliser.newIdentityPreversingUnnormaliser();

	private URL dataUrl;
    private URL indexUrl;
	private ChipsterIndexedFastaSequenceFile picard;

	public IndexedFastaDataSource(URL data, URL index) throws URISyntaxException, IOException {
        super(data, IndexedFastaHandlerThread.class);
        
        this.dataUrl = data;
        this.indexUrl = index;
        
		ChunkDataSource dataSource = new ChunkDataSource(dataUrl, null, null);
		LineDataSource indexSource = new LineDataSource(indexUrl, null);

		picard = new ChipsterIndexedFastaSequenceFile(dataSource, indexSource);
		
		// Create unnormaliser for this naming convention
		this.chromosomeNameUnnormaliser = new ChromosomeNameUnnormaliser(picard.getContigs().get(0));		
    }
	
	public URL getIndex() {
		return indexUrl;
	}

	public URL getUrl() {
		return dataUrl;
	}

	public ChipsterIndexedFastaSequenceFile getPicard() {
		return picard;
	}
	
	public ChromosomeNameUnnormaliser getChromosomeNameUnnormaliser() {
		return chromosomeNameUnnormaliser;
	}
}
