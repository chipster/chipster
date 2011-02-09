package fi.csc.microarray.client.gbrowser.intervals;

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

public class IntervalOperations {

	public static interface PairRule {
		public boolean isPair(BpCoordRegion left, BpCoordRegion right);
	}

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
	 * For each interval in first set, finds intersection with second set.
	 * Returns merged intersections.
	 *  
	 * @param leftIntervals first set 
	 * @param rightIntervals second set
	 * @return intersections
	 */
	public LinkedList<BpCoordRegion> intersect(List<RegionContent> leftIntervals, List<RegionContent> rightIntervals, Long minIntersectionLength, boolean mergeIntersecting) {
		return operate(leftIntervals, rightIntervals, new IntersectingPairRule(minIntersectionLength), EXCLUDE_ORPHAN_POLICY, EXCLUDE_ORPHAN_POLICY, mergeIntersecting ? MERGE_PAIR_POLICY : INTERSECT_PAIR_POLICY, true);
	}


	public LinkedList<BpCoordRegion> subtract(List<RegionContent> leftIntervals, List<RegionContent> rightIntervals, Long minIntersectionLength) {
		return operate(leftIntervals, rightIntervals, new IntersectingPairRule(minIntersectionLength), INCLUDE_ORPHAN_POLICY, EXCLUDE_ORPHAN_POLICY, EXCLUDE_PAIR_POLICY, true);
	}

	public LinkedList<BpCoordRegion> merge(List<RegionContent> leftIntervals, List<RegionContent> rightIntervals, Long minIntersectionLength) {
		return operate(leftIntervals, rightIntervals, new IntersectingPairRule(minIntersectionLength), INCLUDE_ORPHAN_POLICY, INCLUDE_ORPHAN_POLICY, MERGE_PAIR_POLICY, true);
	}

	public void print(Iterable<BpCoordRegion> intervals, OutputStream outputStream) {
		PrintWriter out = new PrintWriter(outputStream);
		for (BpCoordRegion interval : intervals) {
			out.println(interval.toString(true));
		}
		out.flush();
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

	public List<RegionContent> loadFile(File input1) throws FileNotFoundException, IOException {
		ChunkDataSource dataSource = new ChunkDataSource(input1, new BEDParser());
		byte[] fileChunk = dataSource.readAll();
		List<ColumnType> columns = Arrays.asList(new ColumnType[] { ColumnType.CHROMOSOME, ColumnType.BP_START, ColumnType.BP_END});
		return dataSource.getFileParser().getAll(new Chunk(new String(fileChunk)), columns);
	}
	
	public LinkedList<BpCoordRegion> operate(List<RegionContent> leftIntervals, List<RegionContent> rightIntervals, PairRule pairRule, OrphanPolicy leftOrphanPolicy, OrphanPolicy rightOrphanPolicy, PairPolicy pairPolicy, boolean mergeContinous) {
		
		// Initialise collectors
		LinkedList<BpCoordRegion> result = new LinkedList<BpCoordRegion>();
		HashSet<RegionContent> leftPaired = new HashSet<RegionContent>();
		HashSet<RegionContent> rightPaired = new HashSet<RegionContent>();
		
		// Find pairs
		for (RegionContent leftInterval : leftIntervals) {
			for (RegionContent rightInterval : rightIntervals) {
				if (pairRule.isPair(leftInterval.region, rightInterval.region)) {
					leftPaired.add(leftInterval);
					rightPaired.add(rightInterval);
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


}
