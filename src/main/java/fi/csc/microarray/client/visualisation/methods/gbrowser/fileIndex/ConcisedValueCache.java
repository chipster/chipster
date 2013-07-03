package fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex;

import java.util.LinkedHashSet;
import java.util.SortedMap;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;

/**
 * Cache for quickly fetching summary values. By using the cache, we do not have to sample through the file
 * all the time. Works with FIFO principle.  
 * 
 * @author Aleksi Kallio
 *
 */
public class ConcisedValueCache {

	private static final int CACHE_MAX_SIZE = 1000;

	public static class Counts {

		public Integer forwardCount;
		public Integer reverseCount;

		public Counts(int forwardCount, int reverseCount) {
			this.forwardCount = forwardCount;
			this.reverseCount = reverseCount ; 
		}
	}

	/**
	 * Maps location to count of reads on a sample. Same sample size must be used
	 * for all indexed items.
	 */
	private TreeMap<BpCoord, Counts> coordinateOrder = new TreeMap<BpCoord, Counts>();
	
	/**
	 * Linked hash set for keeping track of storage order. Storage works
	 * in FIFO principle, using the linked list. Hashing is needed for 
	 * quickly retrieving intermediate values, because values are touched 
	 * when they are returned. 
	 */
	private LinkedHashSet<BpCoord> storageOrder = new LinkedHashSet<BpCoord>();

	public SortedMap<BpCoord, Counts> subMap(BpCoord from, BpCoord to) {
		
		// Use tree to find the submap
		SortedMap<BpCoord, Counts> subMap = coordinateOrder.subMap(from, to);
		
		// Touch returned values
		for (BpCoord coord : subMap.keySet()) {
			storageOrder.remove(coord);
			storageOrder.add(coord);
		}

		return subMap;
	}

	public void store(BpCoord bpCoord, int countForward, int countReverse) {
		coordinateOrder.put(bpCoord, new Counts(countForward, countReverse));
		storageOrder.add(bpCoord);
		shrink();
	}

	/**
	 * Shrinks cache enough to make it fit to maximum limit. Records are removed in FIFO
	 * order.
	 */
	private void shrink() {
		while (coordinateOrder.size() > CACHE_MAX_SIZE) {
			BpCoord oldest = storageOrder.iterator().next();
			storageOrder.remove(oldest);
			coordinateOrder.remove(oldest);
		}
	}


}