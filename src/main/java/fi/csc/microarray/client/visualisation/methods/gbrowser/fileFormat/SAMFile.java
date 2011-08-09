package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMSequenceRecord;

/**
 * SAM and BAM file abstraction. Uses picard library for actual
 * parsing. Does result splitting and caching.
 * 
 * TODO Check how and when we could close the files.
 * 
 * @author naktinis
 * @see <a href="http://samtools.sourceforge.net/">SAMtools</a>
 * @see fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.SAMFileTest SAMFileTest
 *
 */
public class SAMFile {

	public ChromosomeNameUnnormaliser getChromosomeNameUnnormaliser() {
		return chromosomeNameUnnormaliser;
	}

	public SAMFileReader getReader() {
		return reader;
	}

	private ChromosomeNameUnnormaliser chromosomeNameUnnormaliser = ChromosomeNameUnnormaliser.newIdentityPreversingUnnormaliser();
	private SAMFileReader reader;

	/**
     * @param samFile - SAM or BAM file.
     * @param indexFile - SAM index file (usually with .bai extension).
     */
    public SAMFile(File samFile, File indexFile) {
    	// BAMFileReader emits useless warning to System.err that can't be turned off,
    	// so we direct it to other stream and discard. 
    	PrintStream originalErr = System.err;
    	System.setErr(new PrintStream(new ByteArrayOutputStream()));
        this.reader = new SAMFileReader(samFile, indexFile);

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
    
}
