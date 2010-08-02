package fi.csc.chipster.tools;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;
import java.util.List;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;
import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ElandParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.TsvParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.messaging.JobState;

/**
 * Tool for preprocessing ELAND and BAM data: converts to BAM, sorts and
 * creates an index.
 */
public class NGSPreprocessorTool extends JavaAnalysisJobBase {

	private TsvParser[] parsers = {
			new ElandParser()
	};
	
	@Override
	public String getSADL() {
		
		StringBuffer fileFormats = new StringBuffer();
		for (int i = 0; i < parsers.length; i++) {
			fileFormats.append(parsers[i].getName());
			
			if (i < parsers.length - 1) {
				fileFormats.append(", ");
			}
		}
		
		// TODO more verbose name, name of the second parameter
		return 	"ANALYSIS Utils/Sort (Sort primarily using chromosome and secondarily using start " +
				"location of the feature. File format is used to find columns containing " +
				"chromosome and start location. )" + "\n" +
				
				" INPUT GENERIC input.tsv OUTPUT output.tsv" + "\n" +
				" PARAMETER file.format [" + fileFormats + "] DEFAULT " + parsers[0].getName() + " ()" + "\n";
 	}

	@Override
	protected void execute() { 
		updateState(JobState.RUNNING, "Sorting file");
		
		File inputFile = new File(jobWorkDir, "input.tsv");
		File outputFile = new File(jobWorkDir, "output.tsv");
		
		try {
            elandToBAM(inputFile, outputFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
		
        // TODO sort after creating BAM
//		// get the file format and definitions
//		FileDefinition def = null;
//		for (int i = 0; i < parsers.length; i++) {
//			if (parsers[i].getName().equals(inputMessage.getParameters().get(0))) {
//				def = parsers[i].getFileDefinition();
//			}
//		}		
//
//		// run sorter
//		try {
//			new TsvSorter().sort(inputFile, outputFile, 
//					def.indexOf(ColumnType.CHROMOSOME), def.indexOf(ColumnType.BP_START));
//		} catch (Exception e) {
//			updateState(JobState.FAILED, e.getMessage());
//		}
		
		updateState(JobState.RUNNING, "sort finished");
	}
	
	/**
	 * Convert Eland format to BAM.
	 * 
	 * @param elandFile
	 * @throws IOException 
	 */
	private void elandToBAM(File elandFile, File bamFile) throws IOException {
	    // Input file
	    ChunkDataSource elandData = new ChunkDataSource(elandFile, new ElandParser());
	    
	    // Output file
	    // FIXME: filename, header, bam
	    SAMFileHeader samHeader = new SAMFileHeader();
        SAMFileWriter samWriter = new SAMFileWriterFactory().makeSAMWriter(samHeader, false, bamFile);
	    
	    long filePosition = 0,
	         bytesRead = 1;
	    int increment = 200000;
	    byte[] chunkBytes;
	    while (bytesRead > 0) {
	        // Read some data from ELAND
	        chunkBytes = new byte[increment];
	        bytesRead = elandData.read(filePosition, chunkBytes);
	        String chunkString = new String(chunkBytes);
	        int lastNewLine = chunkString.lastIndexOf("\n");
	        
	        if (lastNewLine == -1) {
	            break;
	        }
	        
	        // Read chunk to a region list
	        Chunk chunk = new Chunk(chunkString.substring(0, lastNewLine));
	        List<RegionContent> contents = elandData.getFileParser().getAll(chunk,
	                Arrays.asList(new ColumnType[] {
	                    ColumnType.ID, ColumnType.CHROMOSOME,
	                    ColumnType.BP_START, ColumnType.STRAND,
	                    ColumnType.SEQUENCE
	                }));
	        
	        // Create SAM records for each alignment
	        for (RegionContent content : contents) {
	            SAMRecord samRecord = new SAMRecord(samHeader);
	            samRecord.setReadName((String)content.values.get(ColumnType.ID));
	            samRecord.setAlignmentStart(((Long)content.values.get(ColumnType.BP_START)).intValue());
	            samRecord.setReadBases(((String)content.values.get(ColumnType.SEQUENCE)).getBytes());
	            samRecord.setReferenceName(((Chromosome)content.values.get(ColumnType.CHROMOSOME)).toString());
	            samRecord.setReadNegativeStrandFlag((Strand)content.values.get(ColumnType.STRAND) == Strand.REVERSED);
	            samRecord.setCigarString(((String)content.values.get(ColumnType.SEQUENCE)).length() + "M");
	            samWriter.addAlignment(samRecord);
            }
	        
	        // Increase the file position
	        filePosition += lastNewLine + 1;
	    }
	    
	    samWriter.close();
	    
	    bamFile.getAbsolutePath();
	}
}
