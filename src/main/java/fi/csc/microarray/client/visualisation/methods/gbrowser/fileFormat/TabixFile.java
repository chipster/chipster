package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.TabixReader;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ConcisedValueCache.Counts;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Tabix file abstraction. Uses tabix-software for actual
 * parsing. This is just a facade.
 * 
 * @author naktinis, akallio, klemela
 *
 */
public class TabixFile {
	
	private ConcisedValueCache cache = new ConcisedValueCache();
	
	private File tabixFile;

    public TabixFile(File tabixFile) {
        this.tabixFile = tabixFile;
    }
    
    /**
     * Find reads in a given range.
     * @throws IOException 
     * 
     */
    public List<RegionContent> getReads(AreaRequest request) throws IOException {

        return TabixReader.query(tabixFile.getAbsolutePath(), "s" /*request.start.chr.toString()*/,
                request.start.bp.toString(), request.end.bp.toString());
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
     * @throws IOException 
     */
    public List<RegionContent> getConciseReads(AreaRequest request) throws IOException {

    	List<RegionContent> responseList = new LinkedList<RegionContent>();
        
        // How many times file is read
        int SAMPLING_GRANULARITY = 40;
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

        	// Use cached content, if exists for this region
        	Collection<Counts> indexedValues = cache.subMap(from, to).values();
        	if (!indexedValues.isEmpty()) {
        		
        		cacheHits++;
        		
        		// sum all
        		for (Counts value : indexedValues) {
        			countForward += value.forwardCount;
        		}
        		
        		// divide to get mean
        		countForward /= indexedValues.size();
        		
        	} else {
        		
        		cacheMisses++;
        		        		
                Collection<RegionContent> queryResult = TabixReader.query(tabixFile.getAbsolutePath(), "s" /*request.start.chr.toString()*/,
                		"" + (stepMiddlepoint - SAMPLE_SIZE/2), "" + (stepMiddlepoint + SAMPLE_SIZE/2));

   				countForward += queryResult.size();

            	// Store value in cache
            	cache.store(new BpCoord((long)stepMiddlepoint, request.start.chr), countForward, 0);
        	}

        	// Create two approximated response objects: one for each strand
        	BpCoordRegion recordRegion =
        		new BpCoordRegion(pos, pos + step, request.start.chr);

        	// Forward
        	HashMap<ColumnType, Object> values = new HashMap<ColumnType, Object>();
        	values.put(ColumnType.VALUE, (float)countForward);
        	values.put(ColumnType.STRAND, Strand.FORWARD);
        	responseList.add(new RegionContent(recordRegion, values));
        }	
        
//        System.out.println("Cache hits: " + cacheHits + ", misses: " + cacheMisses);
        return responseList;
    }
}
