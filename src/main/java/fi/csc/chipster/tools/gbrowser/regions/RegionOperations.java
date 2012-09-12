package fi.csc.chipster.tools.gbrowser.regions;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.net.URISyntaxException;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;

import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.ChunkTreeHandlerThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.BEDParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.util.Strings;

/**
 * A generic tool for operating on Chipster formatted genomic region sets.
 * 
 * @author Aleksi Kallio
 *
 */
public class RegionOperations {

	private static final String EMPTY_EXTRA_FIELDS = "";

	public static void main(String[] args) throws FileNotFoundException, IOException {
		RegionOperations tool = new RegionOperations();
		List<RegionContent> file1 = null;
		List<RegionContent> file2 = null;
		try {
			file1 = tool.loadFile(new File("test1.bed"));
			file2 = tool.loadFile(new File("test2.bed"));
		} catch (URISyntaxException e) {
			e.printStackTrace();
		}

		tool.print(tool.intersect(file1, file2, 1L, RegionOperations.LEFT_PAIR_POLICY_WITH_AUGMENTATION, false), System.out);
	}
	
	
	/**
	 * Intersects regions on two sets and returns either unions or intersections of intersecting pairs.
	 * 
	 * @param leftRegions first set
	 * @param rightRegions second set
	 * @param minIntersectionLength minimum number of shared bases
	 * @param mergeOrIntersect if true return union, otherwise intersection
	 * @return
	 */
	public LinkedList<RegionContent> intersect(List<RegionContent> leftRegions, List<RegionContent> rightRegions, Long minIntersectionLength, PairPolicy pairPolicy, boolean flatten) {
		return operate(leftRegions, rightRegions, new IntersectingPairRule(minIntersectionLength), EXCLUDE_ORPHAN_POLICY, EXCLUDE_ORPHAN_POLICY, pairPolicy, flatten);
	}

	public LinkedList<RegionContent> subtract(List<RegionContent> leftRegions, List<RegionContent> rightRegions, Long minIntersectionLength) {
		return operate(leftRegions, rightRegions, new IntersectingPairRule(minIntersectionLength), INCLUDE_ORPHAN_POLICY, EXCLUDE_ORPHAN_POLICY, EXCLUDE_PAIR_POLICY, true);
	}

	public LinkedList<RegionContent> merge(List<RegionContent> leftRegions, List<RegionContent> rightRegions, Long minIntersectionLength, boolean flatten) {
		return operate(leftRegions, rightRegions, new IntersectingPairRule(minIntersectionLength), INCLUDE_ORPHAN_POLICY, INCLUDE_ORPHAN_POLICY, MERGE_PAIR_POLICY, flatten);
	}

	public LinkedList<RegionContent> flatten(List<RegionContent> leftRegions) {
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
	public LinkedList<RegionContent> operate(List<RegionContent> leftRegions, List<RegionContent> rightRegions, PairRule pairRule, OrphanPolicy leftOrphanPolicy, OrphanPolicy rightOrphanPolicy, PairPolicy pairPolicy, boolean mergeContinous) {
		
		// Initialise collectors
		LinkedList<RegionContent> result = new LinkedList<RegionContent>();
		HashSet<RegionContent> leftPaired = new HashSet<RegionContent>();
		HashSet<RegionContent> rightPaired = new HashSet<RegionContent>();
		
		// Find pairs from a Cartesian product
		for (RegionContent leftRegion : leftRegions) {
			for (RegionContent rightRegion : rightRegions) {
				if (pairRule.isPair(leftRegion.region, rightRegion.region)) {
					leftPaired.add(leftRegion);
					rightPaired.add(rightRegion);
					
					// Output what pair policy dictates
					pairPolicy.process(leftRegion, rightRegion, result);
				}
			}
		}
		
		// Process left orphans
		for (RegionContent leftRegion : leftRegions) {
			if (!leftPaired.contains(leftRegion)) {
				leftOrphanPolicy.process(leftRegion, result);
			}
		}
		
		// Process right orphans
		for (RegionContent rightRegion : rightRegions) {
			if (!rightPaired.contains(rightRegion)) {
				rightOrphanPolicy.process(rightRegion, result);
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
		public boolean isPair(Region left, Region right);
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

		public boolean isPair(Region left, Region right) {
			if (!left.intersects(right)) {
				return false;
			}
			Region intersection = left.intersect(right);
			return intersection.getLength() >= minLength;
		}
	}
	
	/**
	 * Decides what should be output when a pair is found. 
	 */
	public static interface PairPolicy {
		public void process(RegionContent left, RegionContent right, LinkedList<RegionContent> collector);
	}

	public static PairPolicy ORIGINALS_PAIR_POLICY = new PairPolicy() {
		public void process(RegionContent left, RegionContent right, LinkedList<RegionContent> collector) {
			collector.add(left);
			collector.add(right);
		}
	};

	public static PairPolicy MERGE_PAIR_POLICY = new PairPolicy() {
		public void process(RegionContent left, RegionContent right, LinkedList<RegionContent> collector) {
			collector.add(new RegionContent(left.region.merge(right.region), EMPTY_EXTRA_FIELDS));
		}
	};

	public static PairPolicy LEFT_PAIR_POLICY = new PairPolicy() {
		public void process(RegionContent left, RegionContent right, LinkedList<RegionContent> collector) {
			collector.add(left);
		}
	};

	public static PairPolicy LEFT_PAIR_POLICY_WITH_AUGMENTATION = new PairPolicy() {
		public void process(RegionContent left, RegionContent right, LinkedList<RegionContent> collector) {
			collector.add(augment(left, right));
		}
	};

	public static PairPolicy RIGHT_PAIR_POLICY = new PairPolicy() {
		public void process(RegionContent left, RegionContent right, LinkedList<RegionContent> collector) {
			collector.add(right);
		}
	};

	public static PairPolicy RIGHT_PAIR_POLICY_WITH_AUGMENTATION = new PairPolicy() {
		public void process(RegionContent left, RegionContent right, LinkedList<RegionContent> collector) {
			collector.add(augment(right, left));
		}
	};

	public static PairPolicy INTERSECT_PAIR_POLICY = new PairPolicy() {
		public void process(RegionContent left, RegionContent right, LinkedList<RegionContent> collector) {
			collector.add(new RegionContent(left.region.intersect(right.region), EMPTY_EXTRA_FIELDS));
		}
	};

	public static PairPolicy EXCLUDE_PAIR_POLICY = new PairPolicy() {
		public void process(RegionContent left, RegionContent right, LinkedList<RegionContent> collector) {
			// do nothing
		}
	};


	/**
	 * Decides what is done to non-paired regions.
	 */
	public static interface OrphanPolicy {
		public void process(RegionContent region, LinkedList<RegionContent> collector);
	}

	
	public static OrphanPolicy INCLUDE_ORPHAN_POLICY = new OrphanPolicy() {
		public void process(RegionContent region, LinkedList<RegionContent> collector) {
			collector.add(region);
		}
	};

	public static OrphanPolicy EXCLUDE_ORPHAN_POLICY = new OrphanPolicy() {
		public void process(RegionContent region, LinkedList<RegionContent> collector) {
			// do nothing
		}
	};


	/**
	 * Parses regions from a BED text formatted input file.
	 * 
	 * @param input BED file
	 * @return regions and their extra data
	 * @throws URISyntaxException 
	 */
	public List<RegionContent> loadFile(File input) throws FileNotFoundException, IOException, URISyntaxException {
		ChunkDataSource dataSource = new ChunkDataSource(input.toURI().toURL(), new BEDParser(), ChunkTreeHandlerThread.class);
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
		int headerLength = (int)parser.getHeaderLength(string);
		string = string.substring(headerLength > 0 ? (headerLength + 1) : 0);
		
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


	private LinkedList<RegionContent> mergeContinuous(LinkedList<RegionContent> regions) {
		
		// Sort to bring continuous pieces together
		sort(regions);
		
		// Write out continuous regions
		LinkedList<RegionContent> mergedRegions = new LinkedList<RegionContent>();
		for (int i = 0; i < regions.size(); ) {
			
			// Iterate as long as continuous
			int j = i;
			for (; j < regions.size() - 1; ) {
				if (regions.get(i).region.intersects(regions.get(j + 1).region)) {
					// Should be merged, we can continue to look for continuous stuff
					j++;
					
				} else {
					// Here is a gap
					break;
				}
			}
			
			// Write out
			mergedRegions.add(new RegionContent(new Region(regions.get(i).region.start, regions.get(j).region.end), EMPTY_EXTRA_FIELDS));
			
			// Jump to region after the previously written one
			i = j+1;
		}
		
		return mergedRegions;
	}

	/**
	 * Prints out regions and their extra fields in BED text format (without track header row).
	 */
	public void print(List<RegionContent> regionContents, OutputStream outputStream) {
		PrintWriter out = new PrintWriter(outputStream);
		for (RegionContent regionContent : regionContents) {
			out.println(regionContent.toString());
		}
		out.flush();
	}

	/**
	 * Prints out regions and their extra fields in TSV text format. Output has
	 * header row and row name column. Uses extra value ID as row names. 
	 * Row names must be unique, so repeating ID's are ignored.
	 * 
	 */
	public void printTSV(List<RegionContent> regionContents, OutputStream outputStream) {

		// Print header (row name column is nameless)
		PrintWriter out = new PrintWriter(outputStream);
		out.println("chromosome\tstart\tend\t" + Strings.delimit(regionContents.get(0).values.keySet(), "\t").toLowerCase());
		
		// Print row names and values
		HashSet<String> rowNames = new HashSet<String>();
		for (RegionContent regionContent : regionContents) {
			String rowName = regionContent.values.get(ColumnType.ID).toString();
			if (rowNames.contains(rowName)) {
				continue;
			}
			out.println(rowName + "\t" + regionContent.toString());
			rowNames.add(rowName);
		}
		
		out.flush();
	}

	public void sort(List<RegionContent> rows) {
		Collections.sort(rows);
	}
	
	/**
	 * Augment extra fields of primary RegionContent with secondary RegionContent.
	 * @return
	 */
	private static RegionContent augment(RegionContent primary, RegionContent secondary) {
		
		RegionContent augmented = new RegionContent(primary.region, "");
		augmented.values.clear();
		
		// Copy from primary, but augmenting from secondary
		for (Entry<ColumnType, Object> entry : primary.values.entrySet()) {
			Object value = entry.getValue();
			if (value == null || "".equals(value.toString().trim())) {
				Object valueInSecondary = secondary.values.get(entry.getKey());
				if (valueInSecondary != null && !"".equals(valueInSecondary.toString().trim())) {
					augmented.values.put(entry.getKey(), valueInSecondary);
					continue;
				}
			}
			augmented.values.put(entry.getKey(), value);
		}
		
		return augmented;
	}

}
