package fi.csc.chipster.tools.gbrowser.regions;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.BEDParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * A generic tool for operating on Chipster formatted genomic region sets.
 * 
 * @author Aleksi Kallio
 *
 */
public class RegionOperations {

	/**
	 * Intersects regions on two sets and returns either unions or intersections of intersecting pairs.
	 * 
	 * @param leftIntervals first set
	 * @param rightIntervals second set
	 * @param minIntersectionLength minimum number of shared bases
	 * @param mergeOrIntersect if true return union, otherwise intersection
	 * @return
	 */
	public LinkedList<BpCoordRegion> intersect(List<RegionContent> leftIntervals, List<RegionContent> rightIntervals, Long minIntersectionLength, boolean mergeOrIntersect) {
		return operate(leftIntervals, rightIntervals, new IntersectingPairRule(minIntersectionLength), EXCLUDE_ORPHAN_POLICY, EXCLUDE_ORPHAN_POLICY, mergeOrIntersect ? MERGE_PAIR_POLICY : INTERSECT_PAIR_POLICY, true);
	}

	public LinkedList<BpCoordRegion> subtract(List<RegionContent> leftIntervals, List<RegionContent> rightIntervals, Long minIntersectionLength) {
		return operate(leftIntervals, rightIntervals, new IntersectingPairRule(minIntersectionLength), INCLUDE_ORPHAN_POLICY, EXCLUDE_ORPHAN_POLICY, EXCLUDE_PAIR_POLICY, true);
	}

	public LinkedList<BpCoordRegion> merge(List<RegionContent> leftIntervals, List<RegionContent> rightIntervals, Long minIntersectionLength, boolean flatten) {
		return operate(leftIntervals, rightIntervals, new IntersectingPairRule(minIntersectionLength), INCLUDE_ORPHAN_POLICY, INCLUDE_ORPHAN_POLICY, MERGE_PAIR_POLICY, flatten);
	}

	public LinkedList<BpCoordRegion> flatten(List<RegionContent> leftIntervals) {
		return operate(leftIntervals, new LinkedList<RegionContent>(), new IntersectingPairRule(0L), INCLUDE_ORPHAN_POLICY, EXCLUDE_ORPHAN_POLICY, MERGE_PAIR_POLICY, true);
	}

	/**
	 * Generic algorithm for region manipulation. Most of the functionality in this class is based on this method, with different
	 * parameter settings. The algorithm is based on Cartesian product, with configurable pairing rule and additional handling for 
	 * orphan (non-paired) regions.
	 * 
	 * 
	 * @param leftIntervals first set (primary set in some cases)
	 * @param rightIntervals second set
	 * @param pairRule rule for deciding of two regions are a pair
	 * @param leftOrphanPolicy what to do with non-paired regions in first set?
	 * @param rightOrphanPolicy what to do with non-paired regions in second set?
	 * @param pairPolicy how to output a pair of regions
	 * @param mergeContinous if true, continuous pieces of the result are merged in the end
	 * 
	 * @return results produced from paired regions and non-paired regions
	 */
	public LinkedList<BpCoordRegion> operate(List<RegionContent> leftIntervals, List<RegionContent> rightIntervals, PairRule pairRule, OrphanPolicy leftOrphanPolicy, OrphanPolicy rightOrphanPolicy, PairPolicy pairPolicy, boolean mergeContinous) {
		
		// Initialise collectors
		LinkedList<BpCoordRegion> result = new LinkedList<BpCoordRegion>();
		HashSet<RegionContent> leftPaired = new HashSet<RegionContent>();
		HashSet<RegionContent> rightPaired = new HashSet<RegionContent>();
		
		// Find pairs from a Cartesian product
		for (RegionContent leftInterval : leftIntervals) {
			for (RegionContent rightInterval : rightIntervals) {
				if (pairRule.isPair(leftInterval.region, rightInterval.region)) {
					leftPaired.add(leftInterval);
					rightPaired.add(rightInterval);
					
					// Output what pair policy dictates
					pairPolicy.process(leftInterval.region, rightInterval.region, result);
				}
			}
		}
		
		// Process left orphans
		for (RegionContent leftInterval : leftIntervals) {
			if (!leftPaired.contains(leftInterval)) {
				leftOrphanPolicy.process(leftInterval.region, result);
			}
		}
		
		// Process right orphans
		for (RegionContent rightInterval : rightIntervals) {
			if (!rightPaired.contains(rightInterval)) {
				rightOrphanPolicy.process(rightInterval.region, result);
			}
		}
		
		// Do final merge, if needed
		if (mergeContinous) {
			return mergeContinuous(result);
		} else {
			return result;
		}
	}



	/**
	 * Decides if two regions are a pair.
	 */
	public static interface PairRule {
		/**
		 * @param left
		 * @param right
		 * @return true iff left and right are pairs according to this rule
		 */
		public boolean isPair(BpCoordRegion left, BpCoordRegion right);
	}

	/**
	 * Implements pair rule based on intersection. Minimum length can be 
	 * configured so that short overlaps are not accepted.
	 *
	 */
	public static class IntersectingPairRule implements PairRule {

		private Long minLength;

		public IntersectingPairRule(Long minLength) {
			this.minLength = minLength;
		}

		public boolean isPair(BpCoordRegion left, BpCoordRegion right) {
			if (!left.intersects(right)) {
				return false;
			}
			BpCoordRegion intersection = left.intersect(right);
			return intersection.getLength() >= minLength;
		}
	}
	
	/**
	 * Decides what should be output when a pair is found. 
	 */
	public static interface PairPolicy {
		public void process(BpCoordRegion left, BpCoordRegion right, LinkedList<BpCoordRegion> collector);
	}

	public static PairPolicy MERGE_PAIR_POLICY = new PairPolicy() {
		public void process(BpCoordRegion left, BpCoordRegion right, LinkedList<BpCoordRegion> collector) {
			collector.add(left.merge(right));
		}
	};

	public static PairPolicy LEFT_PAIR_POLICY = new PairPolicy() {
		public void process(BpCoordRegion left, BpCoordRegion right, LinkedList<BpCoordRegion> collector) {
			collector.add(left);
		}
	};
	
	public static PairPolicy SUBTRACTED_LEFT_PAIR_POLICY = new PairPolicy() {
		public void process(BpCoordRegion left, BpCoordRegion right, LinkedList<BpCoordRegion> collector) {
			collector.add(left.subtract(right));
		}
	};

	public static PairPolicy RIGHT_PAIR_POLICY = new PairPolicy() {
		public void process(BpCoordRegion left, BpCoordRegion right, LinkedList<BpCoordRegion> collector) {
			collector.add(right);
		}
	};
	
	public static PairPolicy SUBTRACTED_RIGHT_PAIR_POLICY = new PairPolicy() {
		public void process(BpCoordRegion left, BpCoordRegion right, LinkedList<BpCoordRegion> collector) {
			collector.add(right.subtract(left));
		}
	};

	public static PairPolicy INTERSECT_PAIR_POLICY = new PairPolicy() {
		public void process(BpCoordRegion left, BpCoordRegion right, LinkedList<BpCoordRegion> collector) {
			collector.add(left.intersect(right));
		}
	};

	public static PairPolicy EXCLUDE_PAIR_POLICY = new PairPolicy() {
		public void process(BpCoordRegion left, BpCoordRegion right, LinkedList<BpCoordRegion> collector) {
			// do nothing
		}
	};


	/**
	 * Decides what is done to non-paired regions.
	 */
	public static interface OrphanPolicy {
		public void process(BpCoordRegion interval, LinkedList<BpCoordRegion> collector);
	}

	
	public static OrphanPolicy INCLUDE_ORPHAN_POLICY = new OrphanPolicy() {
		public void process(BpCoordRegion interval, LinkedList<BpCoordRegion> collector) {
			collector.add(interval);
		}
	};

	public static OrphanPolicy EXCLUDE_ORPHAN_POLICY = new OrphanPolicy() {
		public void process(BpCoordRegion interval, LinkedList<BpCoordRegion> collector) {
			// do nothing
		}
	};



	/**
	 * Prints out regions in BED-like text format.
	 */
	public void print(Iterable<BpCoordRegion> intervals, OutputStream outputStream) {
		PrintWriter out = new PrintWriter(outputStream);
		for (BpCoordRegion interval : intervals) {
			out.println(interval.toString(true));
		}
		out.flush();
	}

	/**
	 * Parses BED-like region text format from an input file.
	 * 
	 * @param input BED-like file
	 * @return regions and their extra data
	 */
	public List<RegionContent> loadFile(File input) throws FileNotFoundException, IOException {
		ChunkDataSource dataSource = new ChunkDataSource(input, new BEDParser());
		byte[] fileChunk = dataSource.readAll();
		List<ColumnType> columns = Arrays.asList(new ColumnType[] { ColumnType.CHROMOSOME, ColumnType.BP_START, ColumnType.BP_END});
		return dataSource.getFileParser().getAll(new Chunk(new String(fileChunk)), columns);
	}
	

	private LinkedList<BpCoordRegion> mergeContinuous(LinkedList<BpCoordRegion> intervals) {
		
		// Sort to bring continuous pieces together
		Collections.sort(intervals);
		
		// Write out continuous intervals
		LinkedList<BpCoordRegion> mergedIntervals = new LinkedList<BpCoordRegion>();
		for (int i = 0; i < intervals.size(); ) {
			
			// Iterate as long as continuous
			int j = i;
			for (; j < intervals.size() - 1; ) {
				if (intervals.get(i).intersects(intervals.get(j + 1))) {
					// should be merged, we can continue to look for continuous stuff
					j++;
					
				} else {
					// here is a gap
					break;
				}
			}
			
			// Write out
			mergedIntervals.add(new BpCoordRegion(intervals.get(i).start, intervals.get(j).end));
			
			// Jump to interval after the previously written one
			i = j+1;
		}
		
		return mergedIntervals;
	}
}
