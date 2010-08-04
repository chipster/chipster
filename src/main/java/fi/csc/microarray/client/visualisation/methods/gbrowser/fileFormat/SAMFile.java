package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.File;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import net.sf.samtools.BAMFileIndexWriter;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * SAM and BAM file abstraction. Uses picard library for actual
 * parsing. This is just a facade.
 * 
 * @author naktinis
 * @see <a href="http://samtools.sourceforge.net/">SAMtools</a>
 * @see fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.SAMFileTest SAMFileTest
 *
 */
public class SAMFile {
    
    private SAMFileReader reader;
    
    /**
     * FIXME: deprecated, but used in GenomeBrowserStarter.
     * 
     * @param samFile
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
    
    /**
     * @param samFile - SAM or BAM file.
     * @param indexFile - SAM index file (usually with .bai extension).
     */
    public SAMFile(File samFile, File indexFile) {
        reader = new SAMFileReader(samFile, indexFile);
    }
    
    /**
     * Find reads in a given range.
     * 
     * <p>
     * TODO add cigar to the list of returned values
     * <p>
     * TODO add pair information to the list of returned values
     * 
     * @param request
     * @return
     */
    public List<RegionContent> getReads(AreaRequest request) {
        List<RegionContent> responseList = new LinkedList<RegionContent>();
        
        // Read the given region
        CloseableIterator<SAMRecord> iterator =
                reader.query(request.start.chr.toString(),
                request.start.bp.intValue(), request.end.bp.intValue(), false);
        for (Iterator<SAMRecord> i = iterator; i.hasNext();) {
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
            
            // TODO
            
            responseList.add(new RegionContent(recordRegion, values));
        }

        iterator.close();
        return responseList;
    }
    
    /**
     * Return approximation of reads in a given range.
     * <p>
     * Works by dividing requested area into several equally-sized
     * regions (samples) and then for each such region asking only
     * a part of its data (given by SAMPLE_FRACTION) in the beginning
     * of the sample.
     * <p>
     * BAM File requests are still quite expensive.
     * 
     * @param request
     * @return
     */
    public List<RegionContent> getConciseReads(AreaRequest request) {
        List<RegionContent> responseList = new LinkedList<RegionContent>();
        
        // How many times file is read
        int SAMPLE_GRANULARITY = 40;
        
        // How many basepairs are read each time (1/2 means half of sample size)
        float SAMPLE_FRACTION = 1 / 40f;
        
        int step = request.getLength().intValue() / SAMPLE_GRANULARITY;
        int sampleSize = (int)(step * SAMPLE_FRACTION);
        
        for (long pos = request.start.bp; pos < request.end.bp; pos += step) {
            CloseableIterator<SAMRecord> iterator =
                reader.query(request.start.chr.toString(),
                (int)pos, (int)pos + sampleSize, false);
            
            // Count reads in this sample area
            int countReverse = 0;
            int countForward = 0;
            for (Iterator<SAMRecord> i = iterator; i.hasNext();) {
                SAMRecord record = i.next();
                if (record.getReadNegativeStrandFlag()) {
                    countReverse++;
                } else {
                    countForward++;
                }
            }
            
            iterator.close();
            
            // Create two approximated response objects: one for each strand
            BpCoordRegion recordRegion =
                new BpCoordRegion(pos, pos + step, request.start.chr);
            
            // Forward
            HashMap<ColumnType, Object> values = new HashMap<ColumnType, Object>();
            values.put(ColumnType.VALUE, (float)countForward);
            values.put(ColumnType.STRAND, Strand.FORWARD);
            responseList.add(new RegionContent(recordRegion, values));
            
            // Reverse
            values = new HashMap<ColumnType, Object>();
            values.put(ColumnType.VALUE, (float)countReverse);
            values.put(ColumnType.STRAND, Strand.REVERSED);
            responseList.add(new RegionContent(recordRegion, values));
        }
        return responseList;
    }
}
