package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

import net.sf.samtools.BAMFileIndexWriter;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;

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
    
    private SAMFileReader reader;
    
    /**
     * @param samFile - SAM or BAM file.
     */
    public SAMFile(File samFile) {
        reader = new SAMFileReader(samFile);
        if (!reader.hasIndex()) {
            // Create an index file if it is missing
            try {
                // FIXME Number of references
                // FIXME Generated indices don't work
                BAMFileIndexWriter bamIndex = new BAMFileIndexWriter(
                        new File(samFile.getAbsolutePath() + ".bai"), 70000);
                bamIndex.writeBinary(true, samFile.length());
                // Reread the file
                reader.close();
                reader = new SAMFileReader(samFile);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
    
    public List<RegionContent> getReads(BpCoordRegion region) {
        List<RegionContent> responseList = new LinkedList<RegionContent>();
        
        // Read the given region
        CloseableIterator<SAMRecord> iterator =
                reader.query("chr" + region.start.chr.toString(),
                region.start.bp.intValue(), region.end.bp.intValue(), false);
        for(Iterator<SAMRecord> i = iterator; i.hasNext(); ) {
            SAMRecord record = i.next();
            // Region for this read
            BpCoordRegion recordRegion =
                new BpCoordRegion((long)record.getAlignmentStart(),
                        (long)record.getAlignmentEnd(),
                        region.start.chr);
            // Values for this read
            HashMap<ColumnType, Object> values = new HashMap<ColumnType, Object>();
            
            // TODO Deal with "=" and "N" in read string
            values.put(ColumnType.SEQUENCE, record.getReadString());
            values.put(ColumnType.STRAND,
                    record.getReadNegativeStrandFlag() ?
                    Strand.REVERSED : Strand.FORWARD);
            
            responseList.add(new RegionContent(recordRegion, values));
        }

        return responseList;
    }
}
