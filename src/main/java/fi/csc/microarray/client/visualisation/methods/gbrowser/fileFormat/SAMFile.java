package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.File;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import net.sf.samtools.CigarElement;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.util.CloseableIterator;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ConcisedValueCache.Counts;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Cigar;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.CigarItem;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * SAM and BAM file abstraction. Uses picard library for actual
 * parsing. This is just a facade.
 * 
 * TODO Check how and when we could close the files.
 * 
 * @author naktinis
 * @see <a href="http://samtools.sourceforge.net/">SAMtools</a>
 * @see fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.SAMFileTest SAMFileTest
 *
 */
public class SAMFile {
    
	private ConcisedValueCache cache = new ConcisedValueCache();
	public SAMFileReader reader;

	/**
     * @param samFile - SAM or BAM file.
     * @param indexFile - SAM index file (usually with .bai extension).
     */
    public SAMFile(File samFile, File indexFile) {
        this.reader = new SAMFileReader(samFile, indexFile);
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
                this.reader.query(request.start.chr.toString(),
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
            
            if (request.requestedContents.contains(ColumnType.STRAND)) {
            	values.put(ColumnType.STRAND,
            			record.getReadNegativeStrandFlag() ?
            					Strand.REVERSED : Strand.FORWARD);
            	
            }
            
            
            
            if (request.requestedContents.contains(ColumnType.QUALITY)) {
            	
            	/*Now string because of equality problem described below, should be some nice internal
            	 * object type in the future
            	 */
            	
            	values.put(ColumnType.QUALITY, record.getBaseQualityString());            	
            }
            
            if (request.requestedContents.contains(ColumnType.CIGAR)) {      
            	
            	
            	Cigar cigar = new Cigar();
            	
            	for (CigarElement picardElement : record.getCigar().getCigarElements()) {
            		cigar.addElement(new CigarItem(picardElement)); 
            	}
            	            	
            	values.put(ColumnType.CIGAR, cigar);
            }
            
            // TODO Deal with "=" and "N" in read string
            if (request.requestedContents.contains(ColumnType.SEQUENCE)) {            
            	                       	
            	String seq = record.getReadString();
            	
            	//FIXME trying fix changes in reads to fit with the reference sequence, but
            	//this doesn't work yet
//            	StringBuffer buf = new StringBuffer();
//            	
//            	int seqCounter = 0;
//            	
//            	for (CigarElement element : record.getCigar().getCigarElements()) {
//            		if (element.getOperator().consumesReferenceBases()) {
//            			for (int j = 0; j < element.getLength(); j++) {
//            				if (seqCounter + j < seq.length()) {
//            					buf.append(seq.charAt(seqCounter + j));
//            				}
//            				seqCounter++;
//            			}
//            		} else {
//            			for (int j = 0; j < element.getLength(); j++) {
//            				buf.append(" ");
//            			}
//            		}
//            	}
            	
            	values.put(ColumnType.SEQUENCE, seq);
            	//values.put(ColumnType.SEQUENCE, buf.toString());
            }
            
            // TODO Add pair data to values
            
            /* NOTE! RegionContents created from the same read are has to be equal in methods 
             * equals, hash and compareTo. Primary types should be ok,
             * but objects (including tables) has to be handled in those methods separately. 
             * Otherwise tracks keep adding the same reads to their read sets again and again. 
             */
            
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
        int SAMPLING_GRANULARITY = 100;
        int step = request.getLength().intValue() / SAMPLING_GRANULARITY;
        int SAMPLE_SIZE = 100; // FIXME: issue, can be bigger then step size 
        
        int cacheHits = 0;
        int cacheMisses = 0;
        
        for (long pos = request.start.bp; pos < request.end.bp; pos += step) {
        	
        	BpCoord from = new BpCoord(pos, request.start.chr);
        	BpCoord to = new BpCoord(pos + step, request.start.chr);
    		int stepMiddlepoint = (int)pos + step/2;

        	// Count number of reads in a sample from this area
    		int countForward = 0;
    		int countReverse = 0;

        	// Use cached content, if exists for this region
        	Collection<Counts> indexedValues = cache.subMap(from, to).values();
        	if (!indexedValues.isEmpty()) {
        		
        		cacheHits++;
        		
        		// sum all
        		for (Counts value : indexedValues) {
        			countForward += value.forwardCount;
        			countReverse += value.reverseCount;
        		}
        		
        		// divide to get mean
        		countForward /= indexedValues.size();
        		countReverse /= indexedValues.size();
        		
        	} else {
        		
        		cacheMisses++;
        		
        		// Fetch new content by taking sample from the middle of this area
        		int start = stepMiddlepoint - SAMPLE_SIZE/2;
        		int end = stepMiddlepoint + SAMPLE_SIZE/2;
        		CloseableIterator<SAMRecord> iterator =
        			this.reader.query(request.start.chr.toString(),
        					start, end, false);

        		// Count reads in this sample area
        		for (Iterator<SAMRecord> i = iterator; i.hasNext();) {
        			SAMRecord record = i.next();
        			
        			// Accept only records that start in this area (very rough approximation for spliced reads)
        			if (record.getAlignmentStart() >= start && record.getAlignmentEnd() <= end) {
        				if (record.getReadNegativeStrandFlag()) {
        					countReverse++;
        				} else {
        					countForward++;
        				}
        			}
        		}

        		iterator.close();
        		
            	// Store value in cache
            	cache.store(new BpCoord((long)stepMiddlepoint, request.start.chr), countForward, countReverse);
        	}

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
        
//        System.out.println("Cache hits: " + cacheHits + ", misses: " + cacheMisses);
        return responseList;
    }
}
