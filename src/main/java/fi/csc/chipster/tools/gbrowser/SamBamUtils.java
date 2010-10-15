package fi.csc.chipster.tools.gbrowser;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import net.sf.picard.io.IoUtil;
import net.sf.picard.sam.BuildBamIndex;
import net.sf.picard.sam.ViewSam;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;

public class SamBamUtils {
	
	private static final String CHROMOSOME_NAME_PREFIX = "chr";

	public static void convertElandToBam(File elandFile, File bamFile) {
		
	}
	
	public static void sortSamBam(File samBamFile, File sortedBamFile) {
		
		SAMFileReader reader = new SAMFileReader(IoUtil.openFileForReading(samBamFile));
		reader.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.coordinate);
		final SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(reader.getFileHeader(), false, sortedBamFile);
		final Iterator<SAMRecord> iterator = reader.iterator();
		while (iterator.hasNext()) {
			writer.addAlignment(iterator.next());
		}
		reader.close();
		writer.close();
	}

	public static void normaliseBam(File bamFile, File normalisedBamFile) {
		SAMFileReader reader = new SAMFileReader(IoUtil.openFileForReading(bamFile));
		
		SAMSequenceDictionary normalisedDictionary = new SAMSequenceDictionary();
		for (SAMSequenceRecord sequenceRecord : reader.getFileHeader().getSequenceDictionary().getSequences()) {
			String sequenceName = sequenceRecord.getSequenceName();
			
			// Strip prefix
			if (sequenceName.startsWith(CHROMOSOME_NAME_PREFIX)) {
				sequenceName = sequenceName.substring(CHROMOSOME_NAME_PREFIX.length());
			}
			
			normalisedDictionary.addSequence(new SAMSequenceRecord(sequenceName, sequenceRecord.getSequenceLength()));
		}
		reader.getFileHeader().setSequenceDictionary(normalisedDictionary);

		SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(reader.getFileHeader(), true, normalisedBamFile);

		Iterator<SAMRecord> iterator = reader.iterator();
		while (iterator.hasNext()) {
			writer.addAlignment(iterator.next());
		}
		reader.close();
		writer.close();
	}
	
	public static void indexBam(File bamFile, File baiFile) {
		BuildBamIndex.createIndex(new SAMFileReader(IoUtil.openFileForReading(bamFile)), baiFile); 
	}
	
	public static void preprocessSamBam_with_normalise(File samBamFile, File preprocessedBamFile, File baiFile) throws IOException {
		
		// Sort
		File sortedTempBamFile = File.createTempFile("sorted", "bam");
		sortSamBam(samBamFile, sortedTempBamFile);
		
		// Normalise
		normaliseBam(sortedTempBamFile, preprocessedBamFile);
		sortedTempBamFile.delete();

		// Index
		indexBam(preprocessedBamFile, baiFile);
	}

	public static void preprocessSamBam(File samBamFile, File preprocessedBamFile, File baiFile) throws IOException {
		
		// Sort
		sortSamBam(samBamFile, preprocessedBamFile);
		
		// Index
		indexBam(preprocessedBamFile, baiFile);
	}
	
	public static List<String> readChromosomeNames(File bamFile) {
		SAMFileReader reader = null; 
		try {
			reader = new SAMFileReader(IoUtil.openFileForReading(bamFile));
			
			LinkedList<String> chromosomes = new LinkedList<String>();
			for (SAMSequenceRecord record : reader.getFileHeader().getSequenceDictionary().getSequences()) {
				chromosomes.add(record.getSequenceName());
			}
			
			return chromosomes;
			
		} finally {
			closeIfPossible(reader);
		}
		
	}

	private static void closeIfPossible(SAMFileReader reader) {
		if (reader != null) {
			try {
				reader.close();
			} catch (Exception e) {
				// Ignore
			}
		}
	}

	public static void main(String[] args) throws IOException {
		File input = new File("/home/akallio/Desktop/test.sam");
		File temp = File.createTempFile("temp", "bam");
		normaliseBam(input, temp);
		ViewSam.main(new String[] {"I=" + temp.getAbsolutePath()});
	}
}
