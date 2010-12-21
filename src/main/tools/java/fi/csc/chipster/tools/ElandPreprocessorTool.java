package fi.csc.chipster.tools;

import fi.csc.microarray.analyser.JobCancelledException;
import fi.csc.microarray.analyser.java.JavaAnalysisJobBase;

/**
 * Tool for preprocessing ELAND read files: converts to BAM, sorts and
 * creates an index.
 * 
 * TODO: similar tool should be created for processing SAM and BAM files.
 * It should simply sort it by coordinates and create an index. If SAM is
 * given, should also convert to BAM. Also, note the naming of chromosomes.
 * Currently SAMFile assumes that they are named "1", "2" etc.
 */
public class ElandPreprocessorTool extends JavaAnalysisJobBase {

	@Override
	public String getSADL() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	protected void execute() throws JobCancelledException {
		// TODO Auto-generated method stub
		
	}
       
//    int MAX_RECORDS_IN_RAM = 100000;
//    int INDEX_NUM_REFERENCES = 1000;
//
//    private TsvParser[] parsers = {
//            new ElandParser()
//    };
//    
//    
//    /**
//     * Do something when a piece of Eland data is read.
//     */
//    interface ElandReaderHandler {
//        public void process(RegionContent content);
//    }
//    
//    /**
//     * Read chromosome information.
//     */
//    class ElandChrHandler implements ElandReaderHandler {
//        
//        public Set<String> chromosomes = new HashSet<String>(); 
//
//        @Override
//        public void process(RegionContent content) {
//            chromosomes.add(((Chromosome)content.values.get(ColumnType.CHROMOSOME)).toString());
//        }
//        
//    }
//    
//    /**
//     * Read records in an eland file and write to a SAM file.
//     */
//    class ElandSAMHandler implements ElandReaderHandler {
//        
//        private SAMFileWriter samWriter;
//        
//        public ElandSAMHandler(SAMFileWriter samWriter) {
//            this.samWriter = samWriter;
//        }
//        
//        @Override
//        public void process(RegionContent content) {
//            SAMRecord samRecord = new SAMRecord(samWriter.getFileHeader());
//            samRecord.setReadName((String)content.values.get(ColumnType.ID));
//            samRecord.setAlignmentStart(((Long)content.values.get(ColumnType.BP_START)).intValue());
//            samRecord.setReadBases(((String)content.values.get(ColumnType.SEQUENCE)).getBytes());
//            samRecord.setReferenceName(((Chromosome)content.values.get(ColumnType.CHROMOSOME)).toString());
//            samRecord.setReadNegativeStrandFlag((Strand)content.values.get(ColumnType.STRAND) == Strand.REVERSED);
//            samRecord.setCigarString(((String)content.values.get(ColumnType.SEQUENCE)).length() + "M");
//            samWriter.addAlignment(samRecord);
//        }
//        
//    }
//    
//    /**
//     * Read Eland file and process it.
//     */
//    class ElandReader {
//        
//        private ChunkDataSource elandData;
//        private ElandReaderHandler handle;
//        
//        public ElandReader(ChunkDataSource elandData, ElandReaderHandler handle) {
//            this.elandData = elandData;
//            this.handle = handle;
//        }
//        
//        public void read() throws IOException {
//            long filePosition = 0;
//            long bytesRead = 1;
//            int increment = 200000;
//            byte[] chunkBytes;
//            while (bytesRead > 0) {
//                // Read some data from ELAND
//                chunkBytes = new byte[increment];
//                bytesRead = elandData.read(filePosition, chunkBytes);
//                String chunkString = new String(chunkBytes);
//                int lastNewLine = chunkString.lastIndexOf("\n");
//                
//                if (lastNewLine == -1) {
//                    break;
//                }
//                
//                // Read chunk to a region list
//                Chunk chunk = new Chunk(chunkString.substring(0, lastNewLine));
//                List<RegionContent> contents = elandData.getFileParser().getAll(chunk,
//                        Arrays.asList(new ColumnType[] {
//                            ColumnType.ID, ColumnType.CHROMOSOME,
//                            ColumnType.BP_START, ColumnType.STRAND,
//                            ColumnType.SEQUENCE
//                        }));
//                
//                // Somehow process each region content
//                for (RegionContent content : contents) {
//                    handle.process(content);
//                }
//                
//                // Increase the file position
//                filePosition += lastNewLine + 1;
//            }
//        }
//    }
//    
//    @Override
//    public String getSADL() {
//        
//        StringBuffer fileFormats = new StringBuffer();
//        for (int i = 0; i < parsers.length; i++) {
//            fileFormats.append(parsers[i].getName());
//            
//            if (i < parsers.length - 1) {
//                fileFormats.append(", ");
//            }
//        }
//        
//        // TODO more verbose name, name of the second parameter
//        return  "TOOL Utils / eland_preprocessor : \"ELAND preprocessor\"" + 
//                "(Convert ELAND file to BAM, sort it and create an index.)"  + "\n" +
//                "INPUT input.tsv TYPE GENERIC" + "\n" +
//                "OUTPUT output.bam" + "\n" +
//                "OUTPUT output.bai" + "\n";
//     }
//
//    @Override
//    protected void execute() { 
//        updateState(JobState.RUNNING, "Sorting file");
//        
//        File inputFile = new File(jobWorkDir, "input.tsv");
//        File outputFile = new File(jobWorkDir, "output.bam");
//        File indexFile = new File(jobWorkDir, "output.bai");
//        
//        try {
//            elandToBAM(inputFile, outputFile, indexFile);
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
//        
//        updateState(JobState.RUNNING, "sort finished");
//    }
//    
//    /**
//     * Convert Eland format to BAM.
//     * 
//     * @param elandFile
//     * @throws Exception 
//     */
//    private void elandToBAM(File elandFile, File bamFile, File indexFile) throws Exception {
//        // Input file
//        ChunkDataSource elandData = new ChunkDataSource(elandFile, new ElandParser());
//        
//        // Collect chromosomes and later add them to the header
//        ElandChrHandler chrHandler = new ElandChrHandler();
//        new ElandReader(elandData, chrHandler).read();
//        Set<String> chromosomes = chrHandler.chromosomes;
//            
//        // Create header for output file
//        SAMFileHeader samHeader = new SAMFileHeader();
//        samHeader.setSortOrder(SortOrder.coordinate);
//        for (String chr : chromosomes) {
//            samHeader.addSequence(new
//                    SAMSequenceRecord(chr,
//                    (int)GenomeBrowser.CHROMOSOME_SIZES[Integer.parseInt(chr) - 1]));
//        }        
//           
//        // Create SAM writer
//        SAMFileWriterImpl.setDefaultMaxRecordsInRam(MAX_RECORDS_IN_RAM);
//        SAMFileWriter samWriter = new SAMFileWriterFactory().
//                makeBAMWriter(samHeader, false, bamFile);
//        
//        // Write records
//        ElandSAMHandler samHandler = new ElandSAMHandler(samWriter);
//        new ElandReader(elandData, samHandler).read();
//        samWriter.close();
//        
//        // Write index
//        BAMFileIndexWriter bamIndex = new BAMFileIndexWriter(
//                indexFile, INDEX_NUM_REFERENCES);
//        bamIndex.createIndex(bamFile, false, true);
//        bamIndex.writeBinary(true, bamFile.length());      
//        bamIndex.close();
//    }
}
