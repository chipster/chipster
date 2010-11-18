package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.ByteArrayOutputStream;
import java.io.File;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.SortedMap;

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
    	// BAMFileReader emits useless warning to System.err that can't be turned off,
    	// so we direct it to other stream and discard. 
    	PrintStream originalErr = System.err;
    	System.setErr(new PrintStream(new ByteArrayOutputStream()));
        this.reader = new SAMFileReader(samFile, indexFile);
        
        // Restore System.err
        System.setErr(originalErr);
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
            BpCoordRegion recordRegion =
                new BpCoordRegion((long)record.getAlignmentStart(),
                        (long)record.getAlignmentEnd(),
                        request.start.chr);
            
            // Values for this read
            HashMap<ColumnType, Object> values = new HashMap<ColumnType, Object>();
            
            RegionContent read = new RegionContent(recordRegion, values);
            
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
            	
            	Cigar cigar = new Cigar(read);
            	
            	for (CigarElement picardElement : record.getCigar().getCigarElements()) {
            		cigar.addElement(new CigarItem(picardElement)); 
            	}
            	            	
            	values.put(ColumnType.CIGAR, cigar);
            }
            
            // TODO Deal with "=" and "N" in read string
            if (request.requestedContents.contains(ColumnType.SEQUENCE)) {            
            	                       	
            	String seq = record.getReadString();
            	
            	values.put(ColumnType.SEQUENCE, seq);
            }
            
            /* NOTE! RegionContents created from the same read are has to be equal in methods 
             * equals, hash and compareTo. Primary types should be ok,
             * but objects (including tables) has to be handled in those methods separately. 
             * Otherwise tracks keep adding the same reads to their read sets again and again. 
             */
            responseList.add(read);
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
        int SAMPLING_GRANULARITY = 50;
        int step = request.getLength().intValue() / SAMPLING_GRANULARITY;
        
        // Divide visible region into subregions and iterate over them
        for (long pos = request.start.bp; pos < request.end.bp; pos += step) {
        	
        	BpCoord from = new BpCoord(pos, request.start.chr);
        	BpCoord to = new BpCoord(pos + step, request.start.chr);

        	// Use cached content, if exists for this region
        	SortedMap<BpCoord, Counts> indexedValues = cache.subMap(from, to);
        	if (!indexedValues.isEmpty()) {
        		convertCacheHitsToConcisedRegions(request, responseList, step, pos, indexedValues);

        	} else {
        		sampleToGetConcisedRegion(request, responseList, step, pos);
        	}

        }	
        
        return responseList;
    }

	private void sampleToGetConcisedRegion(AreaRequest request, List<RegionContent> responseList, int step, long pos) {

		
        int SAMPLE_SIZE = 100; // FIXME: issue, can be bigger than step size 

		// Fetch new content by taking sample from the middle of this area
		int stepMiddlepoint = (int)pos + step/2;
		int start = stepMiddlepoint - SAMPLE_SIZE/2;
		int end = stepMiddlepoint + SAMPLE_SIZE/2;
		CloseableIterator<SAMRecord> iterator =
			this.reader.query(request.start.chr.toString(),
					start, end, false);

		// Count reads in this sample area
		int countForward = 0;
		int countReverse = 0;
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
		
		addConcisedRegionContents(request, responseList, pos, pos + step, countForward, countReverse);
	}

	private void convertCacheHitsToConcisedRegions(AreaRequest request, List<RegionContent> responseList, int step, long pos, SortedMap<BpCoord, Counts> indexedValues) {
		// Return one result pair for each region covered by one cache hit 
		long startPos = pos;
		int cacheHitsPerRegion = 0;
		for (BpCoord coord : indexedValues.keySet()) {
			
			// Find end: either next cache hit or end of region
			long endPos;
			SortedMap<BpCoord, Counts> tailMap = indexedValues.tailMap(coord);
			if (tailMap.size() > 1) {
				Iterator<BpCoord> iterator = tailMap.keySet().iterator();
				iterator.next(); // read away this cache hit
				BpCoord next = iterator.next(); 
				endPos = (startPos + next.bp) / 2;
				
			} else {
				endPos = (pos + step);
			}
			
			addConcisedRegionContents(request, responseList, startPos, endPos, indexedValues.get(coord).forwardCount, indexedValues.get(coord).reverseCount);
			cacheHitsPerRegion++;
			
			// Move to next region
			startPos += endPos;
		}
//		System.out.println(cacheHitsPerRegion);
	}

	private void addConcisedRegionContents(AreaRequest request, List<RegionContent> responseList, long startPos, long endPos, int countForward, int countReverse) {
		// Create two approximated response objects: one for each strand
		BpCoordRegion recordRegion =
			new BpCoordRegion(startPos, endPos, request.start.chr);

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
}
