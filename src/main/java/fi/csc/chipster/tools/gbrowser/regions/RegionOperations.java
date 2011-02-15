package fi.csc.chipster.tools.gbrowser.regions;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
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
	 * @param leftRegions first set
	 * @param rightRegions second set
	 * @param minIntersectionLength minimum number of shared bases
	 * @param mergeOrIntersect if true return union, otherwise intersection
	 * @return
	 */
	public LinkedList<BpCoordRegion> intersect(List<RegionContent> leftRegions, List<RegionContent> rightRegions, Long minIntersectionLength, boolean mergeOrIntersect) {
		return operate(leftRegions, rightRegions, new IntersectingPairRule(minIntersectionLength), EXCLUDE_ORPHAN_POLICY, EXCLUDE_ORPHAN_POLICY, mergeOrIntersect ? MERGE_PAIR_POLICY : INTERSECT_PAIR_POLICY, true);
	}

	public LinkedList<BpCoordRegion> subtract(List<RegionContent> leftRegions, List<RegionContent> rightRegions, Long minIntersectionLength) {
		return operate(leftRegions, rightRegions, new IntersectingPairRule(minIntersectionLength), INCLUDE_ORPHAN_POLICY, EXCLUDE_ORPHAN_POLICY, EXCLUDE_PAIR_POLICY, true);
	}

	public LinkedList<BpCoordRegion> merge(List<RegionContent> leftRegions, List<RegionContent> rightRegions, Long minIntersectionLength, boolean flatten) {
		return operate(leftRegions, rightRegions, new IntersectingPairRule(minIntersectionLength), INCLUDE_ORPHAN_POLICY, INCLUDE_ORPHAN_POLICY, MERGE_PAIR_POLICY, flatten);
	}

	public LinkedList<BpCoordRegion> flatten(List<RegionContent> leftRegions) {
		return operate(leftRegions, new LinkedList<RegionContent>(), new IntersectingPairRule(0L), INCLUDE_ORPHAN_POLICY, EXCLUDE_ORPHAN_POLICY, MERGE_PAIR_POLICY, true);
	}

	/**
	 * Generic algorithm for region manipulation. Most of the functionality in this class is based on this method, with different
	 * parameter settings. The algorithm is based on Cartesian product, with configurable pairing rule and additional handling for 
	 * orphan (non-paired) regions.
	 * 
	 * 
	 * @param leftRegions first set (primary set in some cases)
	 * @param rightRegions second set
	 * @param pairRule rule for deciding of two regions are a pair
	 * @param leftOrphanPolicy what to do with non-paired regions in first set?
	 * @param rightOrphanPolicy what to do with non-paired regions in second set?
	 * @param pairPolicy how to output a pair of regions
	 * @param mergeContinous if true, continuous pieces of the result are merged in the end
	 * 
	 * @return results produced from paired regions and non-paired regions
	 */
	public LinkedList<BpCoordRegion> operate(List<RegionContent> leftRegions, List<RegionContent> rightRegions, PairRule pairRule, OrphanPolicy leftOrphanPolicy, OrphanPolicy rightOrphanPolicy, PairPolicy pairPolicy, boolean mergeContinous) {
		
		// Initialise collectors
		LinkedList<BpCoordRegion> result = new LinkedList<BpCoordRegion>();
		HashSet<RegionContent> leftPaired = new HashSet<RegionContent>();
		HashSet<RegionContent> rightPaired = new HashSet<RegionContent>();
		
		// Find pairs from a Cartesian product
		for (RegionContent leftRegion : leftRegions) {
			for (RegionContent rightRegion : rightRegions) {
				if (pairRule.isPair(leftRegion.region, rightRegion.region)) {
					leftPaired.add(leftRegion);
					rightPaired.add(rightRegion);
					
					// Output what pair policy dictates
					pairPolicy.process(leftRegion.region, rightRegion.region, result);
				}
			}
		}
		
		// Process left orphans
		for (RegionContent leftRegion : leftRegions) {
			if (!leftPaired.contains(leftRegion)) {
				leftOrphanPolicy.process(leftRegion.region, result);
			}
		}
		
		// Process right orphans
		for (RegionContent rightRegion : rightRegions) {
			if (!rightPaired.contains(rightRegion)) {
				rightOrphanPolicy.process(rightRegion.region, result);
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
		public void process(BpCoordRegion region, LinkedList<BpCoordRegion> collector);
	}

	
	public static OrphanPolicy INCLUDE_ORPHAN_POLICY = new OrphanPolicy() {
		public void process(BpCoordRegion region, LinkedList<BpCoordRegion> collector) {
			collector.add(region);
		}
	};

	public static OrphanPolicy EXCLUDE_ORPHAN_POLICY = new OrphanPolicy() {
		public void process(BpCoordRegion region, LinkedList<BpCoordRegion> collector) {
			// do nothing
		}
	};



	/**
	 * Prints out regions in BED text format (without track header row).
	 */
	public void print(Iterable<BpCoordRegion> regions, OutputStream outputStream) {
		PrintWriter out = new PrintWriter(outputStream);
		for (BpCoordRegion region : regions) {
			out.println(region.toString(true));
		}
		out.flush();
	}

	/**
	 * Parses regions from a BED text formatted input file.
	 * 
	 * @param input BED file
	 * @return regions and their extra data
	 */
	public List<RegionContent> loadFile(File input) throws FileNotFoundException, IOException {
		ChunkDataSource dataSource = new ChunkDataSource(input, new BEDParser());
		byte[] fileChunk = dataSource.readAll();
		return parseString(new String(fileChunk));
	}

	/**
	 * Parses regions from a BED text formatted String.
	 * 
	 * @param string BED string
	 * @return  regions and their extra data
	 */
	public List<RegionContent> parseString(String string) throws FileNotFoundException, IOException {
		
		// Process track name, if exists
		BEDParser parser = new BEDParser();
		string = string.substring((int)parser.getHeaderLength(string) + 1);
		
		// Count fields and create list of what extra types we need
		int fieldCount = string.split("\n")[0].split("\t").length;
		if (fieldCount < 3) {
			throw new IllegalArgumentException("BED must have at least chromosome, start and end fields");
		}
		LinkedList<ColumnType> extraTypes = new LinkedList<ColumnType>();
		for (int i = 3; i < fieldCount; i++) {
			extraTypes.add(BEDParser.completeBedColumns.get(i).content);
		}
		
		// Parse it
		return parser.getAll(new Chunk(string), extraTypes);
	}


	private LinkedList<BpCoordRegion> mergeContinuous(LinkedList<BpCoordRegion> regions) {
		
		// Sort to bring continuous pieces together
		Collections.sort(regions);
		
		// Write out continuous regions
		LinkedList<BpCoordRegion> mergedRegions = new LinkedList<BpCoordRegion>();
		for (int i = 0; i < regions.size(); ) {
			
			// Iterate as long as continuous
			int j = i;
			for (; j < regions.size() - 1; ) {
				if (regions.get(i).intersects(regions.get(j + 1))) {
					// Should be merged, we can continue to look for continuous stuff
					j++;
					
				} else {
					// Here is a gap
					break;
				}
			}
			
			// Write out
			mergedRegions.add(new BpCoordRegion(regions.get(i).start, regions.get(j).end));
			
			// Jump to region after the previously written one
			i = j+1;
		}
		
		return mergedRegions;
	}

	public void printRegions(List<RegionContent> regionContents, OutputStream outputStream) {
		PrintWriter out = new PrintWriter(outputStream);
		for (RegionContent regionContent : regionContents) {
			out.println(regionContent.region.toString(true));
		}
		out.flush();
		
	}
}
