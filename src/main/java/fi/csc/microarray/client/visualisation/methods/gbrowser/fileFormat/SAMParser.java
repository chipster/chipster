package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.File;
import java.util.Collection;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

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
public class SAMParser extends FileParser {
    
    private SAMFileReader reader;
    private SAMFileHeader header;
    
    /**
     * @param samFile - SAM or BAM file.
     */
    public SAMParser(File samFile) {
        reader = new SAMFileReader(samFile);
        header = reader.getFileHeader();
    }

    @Override
    public RegionContent[] concise(String chunk) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public List<RegionContent> getAll(String chunk,
            Collection<ColumnType> requestedContents) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public BpCoordRegion getBpRegion(String chunk) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public long getDefaulChunkLength() {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public String getName() {
        if (reader.isBinary()) {
            return "BAM";
        }
        return "SAM";
    }

}
