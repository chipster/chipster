package fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex;

import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.net.URISyntaxException;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.ChromosomeNameUnnormaliser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.SamBamUtils;
import htsjdk.samtools.SAMFileReader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.util.CloseableIterator;

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
     * @param data
     * @param index
     * @throws URISyntaxException
     * @throws IOException 
     */
    public BamDataSource(DataUrl data, DataUrl index) throws URISyntaxException, IOException {
        super(data);

    	// BAMFileReader emits useless warning to System.err that can't be turned off,
    	// so we direct it to other stream and discard. 
    	PrintStream originalErr = System.err;
    	System.setErr(new PrintStream(new ByteArrayOutputStream()));
    	
    	SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
    	this.reader = SamBamUtils.getSAMReader(data.getUrl(), index.getUrl());

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
