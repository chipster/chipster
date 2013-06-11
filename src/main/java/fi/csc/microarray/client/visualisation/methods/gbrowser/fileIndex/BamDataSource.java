package fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex;

import java.io.ByteArrayOutputStream;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.LinkedList;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.util.CloseableIterator;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.ChromosomeNameUnnormaliser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.SamBamUtils;

/**
 * Data source for indexed BAM compatible files (.bam)
 * 
 * @author Rimvydas Naktinis, Aleksi Kallio
 *
 */
public class BamDataSource extends DataSource {
	
	private ChromosomeNameUnnormaliser chromosomeNameUnnormaliser = ChromosomeNameUnnormaliser.newIdentityPreversingUnnormaliser();
	private SAMFileReader reader;

    /**
     * Generally we would like to have both data and index files,
     * because otherwise we cannot access random locations.
     * 
     * @param samFile
     * @param indexFile
     * @throws FileNotFoundException
     * @throws URISyntaxException
     */
    public BamDataSource(URL samFile, URL indexFile) throws FileNotFoundException, URISyntaxException {
        super(samFile);

    	// BAMFileReader emits useless warning to System.err that can't be turned off,
    	// so we direct it to other stream and discard. 
    	PrintStream originalErr = System.err;
    	System.setErr(new PrintStream(new ByteArrayOutputStream()));
    	
    	this.reader = SamBamUtils.getSAMReader(samFile, indexFile);

    	LinkedList<String> chrList = new LinkedList<>();
    	
    	// Iterate chromosomes to check naming convention
    	for (SAMSequenceRecord sequenceRecord : this.reader.getFileHeader().getSequenceDictionary().getSequences()) {
    		
    		String name = sequenceRecord.getSequenceName();
    		
			chrList.add(name);
    	}
    	
    	// Create unnormaliser for this naming convention
    	// Look only at the first sequence (assume all have the same convention)
    	this.chromosomeNameUnnormaliser = new ChromosomeNameUnnormaliser(chrList);      
        
        // Restore System.err
        System.setErr(originalErr);
    }

	public void close() {
		SamBamUtils.closeIfPossible(reader);
	}	

	public CloseableIterator<SAMRecord> query(Chromosome chr, int start,	int end) {		
		
		String unnormalisedChr = chromosomeNameUnnormaliser.unnormalise(chr);
	
		//The query returns empty collection if there isn't requested chromosome
		SAMRecordIterator iterator = reader.query(unnormalisedChr, start, end, false);

		return iterator; 
	}
}
