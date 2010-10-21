package fi.csc.chipster.tools.gbrowser;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import net.sf.picard.io.IoUtil;
import net.sf.picard.sam.BuildBamIndex;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import fi.csc.microarray.util.IOUtils;

public class SamBamUtils {
	
	private static final String REDUNDANT_CHROMOSOME_NAME_PREFIX = "chr";

	public static void convertElandToSortedBam(File elandFile, File bamFile) throws IOException {

		BufferedReader in = null;
		SAMFileWriter writer = null;
		
		try {
			in = new BufferedReader(new InputStreamReader(new FileInputStream(elandFile)));
			SAMFileHeader header = new SAMFileHeader();
			header.setSortOrder(SAMFileHeader.SortOrder.coordinate);
			writer = new SAMFileWriterFactory().makeBAMWriter(header, false, bamFile);

			for (String line = in.readLine(); line != null; line = in.readLine()) {
				String[] fields = line.split("\t");
				SAMRecord alignment = new SAMRecord(header);
				alignment.setReadName(fields[0]);
				alignment.setReadBases(fields[1].getBytes());
//				String quality =  fields[2];
				alignment.setAlignmentStart(Integer.parseInt(fields[7]));
				alignment.setReadNegativeStrandFlag("R".equals(fields[8]));
				alignment.setReferenceName(fields[6]);
				writer.addAlignment(alignment);
			}
			
		} finally {
			IOUtils.closeIfPossible(in);
			closeIfPossible(writer);
		}

	}
	
	public static void sortSamBam(File samBamFile, File sortedBamFile) {
		
		SAMFileReader reader = new SAMFileReader(IoUtil.openFileForReading(samBamFile));
		SAMFileWriter writer = null;
		try {
			
			reader.getFileHeader().setSortOrder(SAMFileHeader.SortOrder.coordinate);
			writer = new SAMFileWriterFactory().makeBAMWriter(reader.getFileHeader(), false, sortedBamFile);
			Iterator<SAMRecord> iterator = reader.iterator();
			while (iterator.hasNext()) {
				writer.addAlignment(iterator.next());
			}
			
		} finally {
			closeIfPossible(reader);
			closeIfPossible(writer);
		}
	}

	
	public static void normaliseBam(File bamFile, File normalisedBamFile) {

		// Read in a BAM file and its header
		SAMFileReader reader = new SAMFileReader(IoUtil.openFileForReading(bamFile));
		SAMFileWriter writer = null;
		try {
			SAMFileHeader normalisedHeader = reader.getFileHeader();

			// Alter the chromosome names in header's SAMSequenceDictionary
			SAMSequenceDictionary normalisedDictionary = new SAMSequenceDictionary();
			for (SAMSequenceRecord sequenceRecord : normalisedHeader.getSequenceDictionary().getSequences()) {

				// Strip prefix, if exists
				String sequenceName = sequenceRecord.getSequenceName();
				if (sequenceName.startsWith(REDUNDANT_CHROMOSOME_NAME_PREFIX)) {
					sequenceName = sequenceName.substring(REDUNDANT_CHROMOSOME_NAME_PREFIX.length());
				}

				normalisedDictionary.addSequence(new SAMSequenceRecord(sequenceName, sequenceRecord.getSequenceLength()));
			}
			normalisedHeader.setSequenceDictionary(normalisedDictionary);

			// Write new BAM file with normalised chromosome names
			writer = new SAMFileWriterFactory().makeBAMWriter(normalisedHeader, true, normalisedBamFile);
			for (final SAMRecord rec : reader) {
				rec.setHeader(normalisedHeader);
				writer.addAlignment(rec);
			}
			
		} finally {
			closeIfPossible(reader);
			closeIfPossible(writer);
		}
	}

	public static void indexBam(File bamFile, File baiFile) {
		BuildBamIndex.createIndex(new SAMFileReader(IoUtil.openFileForReading(bamFile)), baiFile); 
	}

	public static void preprocessEland(File elandFile, File preprocessedBamFile, File baiFile) throws IOException {
		File sortedTempBamFile = File.createTempFile("converted", "bam");
		
		try {
			// Convert & Sort
			convertElandToSortedBam(elandFile, sortedTempBamFile);

			// Normalise (input must be BAM)
			normaliseBam(sortedTempBamFile, preprocessedBamFile);
			sortedTempBamFile.delete();

			// Index
			indexBam(preprocessedBamFile, baiFile);
			
		} finally {
			sortedTempBamFile.delete();
		}
	}

	public static void preprocessSamBam(File samBamFile, File preprocessedBamFile, File baiFile) throws IOException {
		
		// Sort
		File sortedTempBamFile = File.createTempFile("sorted", "bam");
		sortSamBam(samBamFile, sortedTempBamFile);
		
		// Normalise (input must be BAM)
		normaliseBam(sortedTempBamFile, preprocessedBamFile);
		sortedTempBamFile.delete();

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

	private static void closeIfPossible(SAMFileWriter writer) {
		if (writer != null) {
			try {
				writer.close();
			} catch (Exception e) {
				// Ignore
			}
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

	public static boolean isSamBamExtension(String extension) {
		if (extension == null) {
			return false;
		}
		extension = extension.toLowerCase();
		return "sam".equals(extension) || "bam".equals(extension);
	}

}
