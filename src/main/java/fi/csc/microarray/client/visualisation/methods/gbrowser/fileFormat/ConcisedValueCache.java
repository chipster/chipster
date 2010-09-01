package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import java.util.SortedMap;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;

/**
 * Index for quickly fetching summary values. 
 * 
 * @author Aleksi Kallio
 *
 */
public class ConcisedValueCache {

	public static class Counts {

		Integer forwardCount;
		Integer reverseCount;

		public Counts(int forwardCount, int reverseCount) {
			this.forwardCount = forwardCount;
			this.reverseCount = reverseCount ; 
		}
	}

	/**
	 * Maps location to count of reads on a sample. Same sample size must be used
	 * for all indexed items.
	 */
	public TreeMap<BpCoord, Counts> tree = new TreeMap<BpCoord, Counts>();

	public SortedMap<BpCoord, Counts> subMap(BpCoord from, BpCoord to) {
		return tree.subMap(from, to);
	}

	public void store(BpCoord bpCoord, int countForward, int countReverse) {
		tree.put(bpCoord, new Counts(countForward, countReverse));
	}


}