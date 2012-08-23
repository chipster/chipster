package fi.csc.microarray.client.visualisation.methods.gbrowser;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.net.URISyntaxException;
import java.net.URL;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMSequenceRecord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.SAMHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ChromosomeNameUnnormaliser;

/**
 * Data source for indexed SAM compatible files (.sam/.bam). Because index is required, in practice
 * supports only .bam.
 * 
 * @author Rimvydas Naktinis, Aleksi Kallio
 *
 */
public class SAMDataSource extends DataSource {
	
	private ChromosomeNameUnnormaliser chromosomeNameUnnormaliser = ChromosomeNameUnnormaliser.newIdentityPreversingUnnormaliser();
	private SAMFileReader reader;

    /**
     * Generally we would like to have both data and index files,
     * because otherwise we could not access random locations.
     * 
     * @param file
     * @throws FileNotFoundException
     * @throws URISyntaxException 
     */
    public SAMDataSource(URL samFile, URL indexFile) throws FileNotFoundException, URISyntaxException {
        super(samFile, SAMHandlerThread.class);

    	// BAMFileReader emits useless warning to System.err that can't be turned off,
    	// so we direct it to other stream and discard. 
    	PrintStream originalErr = System.err;
    	System.setErr(new PrintStream(new ByteArrayOutputStream()));
        this.reader = new SAMFileReader(this.file, new File(indexFile.toURI()));

        // Iterate chromosomes to check naming convention
        for (SAMSequenceRecord sequenceRecord : this.reader.getFileHeader().getSequenceDictionary().getSequences()) {

        	// Create unnormaliser for this naming convention
        	this.chromosomeNameUnnormaliser = new ChromosomeNameUnnormaliser(sequenceRecord.getSequenceName());

        	// Look only at the first sequence (assume all have the same convention)
        	break;
        }

        
        // Restore System.err
        System.setErr(originalErr);

    }

	public SAMFileReader getReader() {
		return reader;
	}

	public ChromosomeNameUnnormaliser getChromosomeNameUnnormaliser() {
		return chromosomeNameUnnormaliser;
	}

}
