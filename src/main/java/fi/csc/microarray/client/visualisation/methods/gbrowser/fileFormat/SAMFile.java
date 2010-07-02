package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
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
                // FIXME Slow...
                BAMFileIndexWriter bamIndex = new BAMFileIndexWriter(
                        new File(samFile.getAbsolutePath() + ".bai"), 1000);
                bamIndex.createIndex(samFile, false, true);
                bamIndex.writeBinary(true, samFile.length());
                bamIndex.close();
                // Reread the file
                reader.close();
                reader = new SAMFileReader(samFile);
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
    }
    
    public List<RegionContent> getReads(AreaRequest request) {
        List<RegionContent> responseList = new LinkedList<RegionContent>();
        
        // Read the given region
        CloseableIterator<SAMRecord> iterator =
                reader.query("chr" + request.start.chr.toString(),
                request.start.bp.intValue(), request.end.bp.intValue(), false);
        for(Iterator<SAMRecord> i = iterator; i.hasNext(); ) {
            SAMRecord record = i.next();
            // Region for this read
            // FIXME What happens if it spans across several chromosomes?
            BpCoordRegion recordRegion =
                new BpCoordRegion((long)record.getAlignmentStart(),
                        (long)record.getAlignmentEnd(),
                        request.start.chr);
            // Values for this read
            HashMap<ColumnType, Object> values = new HashMap<ColumnType, Object>();
            
            // TODO Deal with "=" and "N" in read string
            if (request.requestedContents.contains(ColumnType.SEQUENCE)) {
                values.put(ColumnType.SEQUENCE, record.getReadString());
            }
            
            if (request.requestedContents.contains(ColumnType.STRAND)) {
                values.put(ColumnType.STRAND,
                        record.getReadNegativeStrandFlag() ?
                        Strand.REVERSED : Strand.FORWARD);
            }
            
            // TODO Add cigar data to values
            
            responseList.add(new RegionContent(recordRegion, values));
        }

        iterator.close();
        return responseList;
    }
}
