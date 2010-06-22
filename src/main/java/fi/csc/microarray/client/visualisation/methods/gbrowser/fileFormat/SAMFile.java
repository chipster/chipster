package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.File;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;

/**
 * SAM and BAM file abstraction. Uses picard library for actual
 * parsing. This is just an interface.
 * 
 * @author naktinis
 * @see <a href="http://samtools.sourceforge.net/">SAMtools</a>
 * 
 * TODO
 *
 */
public class SAMFile {
    
    public SAMFile(File samFile) {
        SAMFileReader reader = new SAMFileReader(samFile);
        SAMFileHeader header = reader.getFileHeader();
    }
}
