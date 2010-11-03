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
	
	public interface SamBamUtilStateListener {
		public void stateChanged(SamBamUtilState newState);
	}
	
	public class SamBamUtilState {
		private String state;
		private double percentage;
	
		public SamBamUtilState(String state, double percentage) {
			this.state = state;
			this.percentage = percentage;
		}

		public String getState() {
			return this.state;
		}
		
		public double getPercentage() {
			return this.percentage;
		}
	}
	
	
	private static final String CHROMOSOME_NAME_PREFIX = "chr";
	private SamBamUtilStateListener stateListener;
	
	public SamBamUtils() {
	}
	
	public SamBamUtils(SamBamUtilStateListener stateListener) {
		this.stateListener = stateListener;
	}

	private void updateState(String state, double percentage) {
		if (this.stateListener != null) {
			stateListener.stateChanged(new SamBamUtilState(state, percentage));
		}
	}
	
	
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
		
		// Read in a BAM file and its header
		SAMFileReader reader = new SAMFileReader(IoUtil.openFileForReading(bamFile));
		SAMFileHeader normalisedHeader = reader.getFileHeader();
		
		// Alter the chrosome names in header's SAMSequenceDictionary
		SAMSequenceDictionary normalisedDictionary = new SAMSequenceDictionary();
		for (SAMSequenceRecord sequenceRecord : normalisedHeader.getSequenceDictionary().getSequences()) {
			
			// Strip prefix, if exists
			String sequenceName = sequenceRecord.getSequenceName();
			if (sequenceName.startsWith(CHROMOSOME_NAME_PREFIX)) {
				sequenceName = sequenceName.substring(CHROMOSOME_NAME_PREFIX.length());
			}
			
			normalisedDictionary.addSequence(new SAMSequenceRecord(sequenceName, sequenceRecord.getSequenceLength()));
		}
		normalisedHeader.setSequenceDictionary(normalisedDictionary);
		
		// Write new BAM file with normalised chromosome names
		SAMFileWriter writer = new SAMFileWriterFactory().makeBAMWriter(normalisedHeader, true, normalisedBamFile);
		for (final SAMRecord rec : reader) {
			rec.setHeader(normalisedHeader);
			writer.addAlignment(rec);
		}
		
		writer.close();
		reader.close();
	}

	public static void indexBam(File bamFile, File baiFile) {
		BuildBamIndex.createIndex(new SAMFileReader(IoUtil.openFileForReading(bamFile)), baiFile); 
	}
	
	public void preprocessSamBam(File samBamFile, File preprocessedBamFile, File baiFile) throws IOException {
		
		// Sort
		updateState("sorting", 0);
		File sortedTempBamFile = File.createTempFile("sorted", "bam");
		sortSamBam(samBamFile, sortedTempBamFile);
		
		// Normalise (input must be BAM)
		updateState("normalising", 33.3);
		normaliseBam(sortedTempBamFile, preprocessedBamFile);
		sortedTempBamFile.delete();

		// Index
		updateState("indexing", 66);
		indexBam(preprocessedBamFile, baiFile);
		updateState("done", 100);
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
		File input = new File("/home/akallio/Desktop/test.bam");
		File temp = File.createTempFile("temp", "bam");
		normaliseBam(input, temp);
		ViewSam.main(new String[] {"I=" + temp.getAbsolutePath(), "VALIDATION_STRINGENCY=SILENT"});
	}
}
