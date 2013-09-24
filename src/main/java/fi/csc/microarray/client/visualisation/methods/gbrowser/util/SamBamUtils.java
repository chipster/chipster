package fi.csc.microarray.client.visualisation.methods.gbrowser.util;

import java.io.BufferedReader;
import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.net.URISyntaxException;
import java.net.URL;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import net.sf.picard.io.IoUtil;
import net.sf.picard.sam.BuildBamIndex;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceDictionary;
import net.sf.samtools.SAMSequenceRecord;
import net.sf.samtools.seekablestream.SeekableBufferedStream;
import net.sf.samtools.seekablestream.SeekableFileStream;
import net.sf.samtools.seekablestream.SeekableHTTPStream;
import net.sf.samtools.seekablestream.SeekableStream;

import org.broad.tribble.readers.TabixReader;

import fi.csc.microarray.util.IOUtils;

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
	
	
	private SamBamUtilStateListener stateListener;
	private ChromosomeNormaliser chromosomeNormaliser = new ChromosomeNormaliser() {

		public String normaliseChromosome(String chromosomeName) {

			// Leave prefix as it is
			
			// Remove postfix, if present
			String SEPARATOR = ".";
			if (chromosomeName.contains(SEPARATOR)) {
				chromosomeName = chromosomeName.substring(0, chromosomeName.indexOf(SEPARATOR));
			}
			
			return chromosomeName;
		}
	};
	
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
		
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
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

	
	public void normaliseBam(File bamFile, File normalisedBamFile) {

		// Read in a BAM file and its header
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		SAMFileReader reader = new SAMFileReader(IoUtil.openFileForReading(bamFile));
		SAMFileWriter writer = null;
		try {
			SAMFileHeader normalisedHeader = reader.getFileHeader();

			// Alter the chromosome names in header's SAMSequenceDictionary
			SAMSequenceDictionary normalisedDictionary = new SAMSequenceDictionary();
			for (SAMSequenceRecord sequenceRecord : normalisedHeader.getSequenceDictionary().getSequences()) {

				// Normalise chromosome
				String sequenceName = chromosomeNormaliser.normaliseChromosome(sequenceRecord.getSequenceName());
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

	public void indexBam(File bamFile, File baiFile) {
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		BuildBamIndex.createIndex(new SAMFileReader(IoUtil.openFileForReading(bamFile)), baiFile); 
	}

	public void preprocessEland(File elandFile, File preprocessedBamFile, File baiFile) throws IOException {
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

	public static List<String> readChromosomeNames(URL bam, URL index) throws FileNotFoundException, URISyntaxException {
									
			SAMFileReader reader = getSAMReader(bam, index);
			
			LinkedList<String> chromosomes = new LinkedList<String>();
			for (SAMSequenceRecord record : reader.getFileHeader().getSequenceDictionary().getSequences()) {
				chromosomes.add(record.getSequenceName());
			}
			
			closeIfPossible(reader);
			
			return chromosomes;		
	}

	public static SAMFileReader getSAMReader(URL bam, URL index) throws FileNotFoundException, URISyntaxException {
		
		SAMFileReader reader = null; 
		SeekableStream bamStream = null;
		SeekableStream indexStream = null;
		
		try {
			
			if ("file".equals(bam.getProtocol())) {
			
				bamStream = new SeekableFileStream(new File(bam.toURI()));
			} else {				
				bamStream = new SeekableHTTPStream(bam);
			}

			// Picard requires file to have a 'bam' extension, but our hashed url doesn't have it
			SeekableBufferedStream bamBufferedStream = new SeekableBufferedStream(bamStream) {
				@Override
				public String getSource() {
					return super.getSource() + "_fake-source.bam";
				}
			};

			if ("file".equals(index.getProtocol())) {
				
				indexStream = new SeekableFileStream(new File(index.toURI()));
			} else {				
				indexStream = new SeekableHTTPStream(index);
			}
			
			SeekableBufferedStream indexBufferedStream = new SeekableBufferedStream(indexStream);

			return new SAMFileReader(bamBufferedStream, indexBufferedStream, false);

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

	public static void closeIfPossible(SAMFileReader reader) {
		if (reader != null) {
			try {
				reader.close();
			} catch (Exception e) {
				// Ignore
			}
		}
	}
	

	public static void closeIfPossible(TabixReader reader) {
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
	
	public String printSamBam(InputStream samBamStream, int maxRecords) throws IOException {
		SAMFileReader.setDefaultValidationStringency(ValidationStringency.SILENT);
		SAMFileReader in = new SAMFileReader(samBamStream);
		SAMFileHeader header = in.getFileHeader();
		ByteArrayOutputStream buffer = new ByteArrayOutputStream();
		SAMFileWriter out = new SAMFileWriterFactory().makeSAMWriter(header, true, buffer);
		int i = 0;
		try {
			for (final SAMRecord rec : in) {
				if (i > maxRecords) {
					break;
				}
				out.addAlignment(rec);
				i++;
			}
		} finally {
			closeIfPossible(out);
		}

		if (i > maxRecords) {
			buffer.write("SAM/BAM too long for viewing, truncated here!\n".getBytes());
		}
		
		return buffer.toString();
	}
}
