package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.File;
import java.util.Arrays;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;

import org.testng.Assert;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FsfStatus;

public class SAMFileTest {
    private static String path = "src/test/resources/sam/";

    public static void main(String[] args) throws Exception {
        
        // Parse BAM file header using picard library
        File bamFile = new File(path + "variousAttributes.bam");
        SAMFileReader reader = new SAMFileReader(bamFile);
        SAMFileHeader header = reader.getFileHeader();
        Assert.assertEquals(header.getAttribute("VN"), "1.0");
        Assert.assertEquals(header.getSequenceDictionary().getSequence(0).
                getSequenceName(), "chr20");
        
        // Parse BAM file that has index
        bamFile = new File(path + "indexTest.bam");
        reader = new SAMFileReader(bamFile);
        reader.enableIndexCaching(true);
        Assert.assertEquals(reader.getFileHeader().getSequenceDictionary().
                getSequence(0).getSequenceName(), "chrM");
        Assert.assertEquals(reader.query("chrM", 0, 300000, true).next().
                getCigarString(), "51M");
        
        // Parse using our internal abstraction
        bamFile = new File(path + "indexTest.bam");
        File indexFile = new File(path + "indexTest.bam.bai");
        SAMFile sam = new SAMFile(bamFile, indexFile);
        BpCoordRegion region = new BpCoordRegion((long)0, (long)300000,
                new Chromosome("chrM"));
        Assert.assertEquals(
                sam.getReads(new AreaRequest(region,
                        Arrays.asList(new ColumnType[] {ColumnType.SEQUENCE, ColumnType.STRAND}),
                        new FsfStatus())).size(), 23);
    }
}
