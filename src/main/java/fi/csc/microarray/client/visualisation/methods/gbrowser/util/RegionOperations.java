package fi.csc.microarray.client.visualisation.methods.gbrowser.util;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintWriter;
import java.net.URISyntaxException;
import java.util.Collections;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;
import java.util.Map.Entry;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.DataUrl;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.BedLineParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.RandomAccessLineReader;
import fi.csc.microarray.util.Strings;

/**
 * A generic tool for operating on Chipster formatted genomic region sets.
 * 
 * @author Aleksi Kallio
 *
 */
public class RegionOperations {

	private static final String EMPTY_EXTRA_FIELDS = "";

	public static void main(String[] args) throws FileNotFoundException, IOException, GBrowserException {
		RegionOperations tool = new RegionOperations();
		List<Feature> file1 = null;
		List<Feature> file2 = null;
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
	public LinkedList<Feature> intersect(List<Feature> leftRegions, List<Feature> rightRegions, Long minIntersectionLength, PairPolicy pairPolicy, boolean flatten) {
		return operate(leftRegions, rightRegions, new IntersectingPairRule(minIntersectionLength), EXCLUDE_ORPHAN_POLICY, EXCLUDE_ORPHAN_POLICY, pairPolicy, flatten);
	}

	public LinkedList<Feature> subtract(List<Feature> leftRegions, List<Feature> rightRegions, Long minIntersectionLength) {
		return operate(leftRegions, rightRegions, new IntersectingPairRule(minIntersectionLength), INCLUDE_ORPHAN_POLICY, EXCLUDE_ORPHAN_POLICY, EXCLUDE_PAIR_POLICY, true);
	}

	public LinkedList<Feature> merge(List<Feature> leftRegions, List<Feature> rightRegions, Long minIntersectionLength, boolean flatten) {
		return operate(leftRegions, rightRegions, new IntersectingPairRule(minIntersectionLength), INCLUDE_ORPHAN_POLICY, INCLUDE_ORPHAN_POLICY, MERGE_PAIR_POLICY, flatten);
	}

	public LinkedList<Feature> flatten(List<Feature> leftRegions) {
		return operate(leftRegions, new LinkedList<Feature>(), new IntersectingPairRule(0L), INCLUDE_ORPHAN_POLICY, EXCLUDE_ORPHAN_POLICY, MERGE_PAIR_POLICY, true);
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
	public LinkedList<Feature> operate(List<Feature> leftRegions, List<Feature> rightRegions, PairRule pairRule, OrphanPolicy leftOrphanPolicy, OrphanPolicy rightOrphanPolicy, PairPolicy pairPolicy, boolean mergeContinous) {
		
		// Initialise collectors
		LinkedList<Feature> result = new LinkedList<Feature>();
		HashSet<Feature> leftPaired = new HashSet<Feature>();
		HashSet<Feature> rightPaired = new HashSet<Feature>();
		
		// Find pairs from a Cartesian product
		for (Feature leftRegion : leftRegions) {
			for (Feature rightRegion : rightRegions) {
				if (pairRule.isPair(leftRegion.region, rightRegion.region)) {
					leftPaired.add(leftRegion);
					rightPaired.add(rightRegion);
					
					// Output what pair policy dictates
					pairPolicy.process(leftRegion, rightRegion, result);
				}
			}
		}
		
		// Process left orphans
		for (Feature leftRegion : leftRegions) {
			if (!leftPaired.contains(leftRegion)) {
				leftOrphanPolicy.process(leftRegion, result);
			}
		}
		
		// Process right orphans
		for (Feature rightRegion : rightRegions) {
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
		public void process(Feature left, Feature right, LinkedList<Feature> collector);
	}

	public static PairPolicy ORIGINALS_PAIR_POLICY = new PairPolicy() {
		public void process(Feature left, Feature right, LinkedList<Feature> collector) {
			collector.add(left);
			collector.add(right);
		}
	};

	public static PairPolicy MERGE_PAIR_POLICY = new PairPolicy() {
		public void process(Feature left, Feature right, LinkedList<Feature> collector) {
			collector.add(new Feature(left.region.merge(right.region), getEmptyExtraFieldMap()));
		}
	};

	public static PairPolicy LEFT_PAIR_POLICY = new PairPolicy() {
		public void process(Feature left, Feature right, LinkedList<Feature> collector) {
			collector.add(left);
		}
	};

	public static PairPolicy LEFT_PAIR_POLICY_WITH_AUGMENTATION = new PairPolicy() {
		public void process(Feature left, Feature right, LinkedList<Feature> collector) {
			collector.add(augment(left, right));
		}
	};

	public static PairPolicy RIGHT_PAIR_POLICY = new PairPolicy() {
		public void process(Feature left, Feature right, LinkedList<Feature> collector) {
			collector.add(right);
		}
	};

	public static PairPolicy RIGHT_PAIR_POLICY_WITH_AUGMENTATION = new PairPolicy() {
		public void process(Feature left, Feature right, LinkedList<Feature> collector) {
			collector.add(augment(right, left));
		}
	};

	public static PairPolicy INTERSECT_PAIR_POLICY = new PairPolicy() {
		public void process(Feature left, Feature right, LinkedList<Feature> collector) {
			collector.add(new Feature(left.region.intersect(right.region), getEmptyExtraFieldMap()));
		}
	};

	public static PairPolicy EXCLUDE_PAIR_POLICY = new PairPolicy() {
		public void process(Feature left, Feature right, LinkedList<Feature> collector) {
			// do nothing
		}
	};


	/**
	 * Decides what is done to non-paired regions.
	 */
	public static interface OrphanPolicy {
		public void process(Feature region, LinkedList<Feature> collector);
	}

	
	public static OrphanPolicy INCLUDE_ORPHAN_POLICY = new OrphanPolicy() {
		public void process(Feature region, LinkedList<Feature> collector) {
			collector.add(region);
		}
	};

	public static OrphanPolicy EXCLUDE_ORPHAN_POLICY = new OrphanPolicy() {
		public void process(Feature region, LinkedList<Feature> collector) {
			// do nothing
		}
	};


	/**
	 * Parses regions from a BED text formatted input file.
	 * 
	 * @param input BED file
	 * @return regions and their extra data
	 * @throws URISyntaxException 
	 * @throws GBrowserException 
	 */
	public List<Feature> loadFile(File input) throws FileNotFoundException, IOException, URISyntaxException, GBrowserException {
		
		DataUrl dataUrl = new DataUrl(input);
		RandomAccessLineReader lineReader = new RandomAccessLineReader(dataUrl);
		lineReader.setPosition(0);
		
		String line;
		List<Feature> regions = new LinkedList<Feature>();
		
		while ((line = lineReader.readLine()) != null) {	
			
			Feature region = parseString(line);
			
			if (region != null) {
				regions.add(region);
			}
		}
		
		lineReader.close();
		
		return regions;
	}

	/**
	 * Parses region from a BED text formatted String.
	 * 
	 * @param string BED string
	 * @return  regions and their extra data
	 */
	public Feature parseString(String string) throws FileNotFoundException, IOException {
		
		// Process track name, if exists
		BedLineParser parser = new BedLineParser(false);
		parser.setLine(string);
		Region region = parser.getRegion();
		
		if (region == null) {
			//header
			return null;
		}
						
		// Count fields and create list of what extra types we need
		int fieldCount = parser.getColumnCount();
		if (fieldCount < 3) {
			throw new IllegalArgumentException("BED must have at least chromosome, start and end fields");
		}
		
		DataType[] legacyBedColumns = new DataType[] {
				DataType.CHROMOSOME,
				DataType.START,
				DataType.END,
				DataType.ID,
				DataType.VALUE,
				DataType.STRAND,
				DataType.THICK_START,
				DataType.THICK_END,
				DataType.ITEM_RGB,
				DataType.BLOCK_COUNT,
				DataType.BLOCK_SIZES,
				DataType.BLOCK_STARTS
		};
		
		LinkedHashMap<DataType, Object> values = new LinkedHashMap<DataType, Object>();
		
		for (int i = 3; i < fieldCount; i++) {
						
			DataType key = legacyBedColumns[i];
			String value = parser.getString(i);
			values.put(key, value);
		}
		
		return new Feature(region, values);
	}


	private LinkedList<Feature> mergeContinuous(LinkedList<Feature> regions) {
		
		// Sort to bring continuous pieces together
		sort(regions);
		
		// Write out continuous regions
		LinkedList<Feature> mergedRegions = new LinkedList<Feature>();
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
			mergedRegions.add(new Feature(new Region(regions.get(i).region.start, regions.get(j).region.end), getEmptyExtraFieldMap()));
			
			// Jump to region after the previously written one
			i = j+1;
		}
		
		return mergedRegions;
	}

	/**
	 * Prints out regions and their extra fields in BED text format (without track header row).
	 */
	public void print(List<Feature> regionContents, OutputStream outputStream) {
		PrintWriter out = new PrintWriter(outputStream);
		for (Feature regionContent : regionContents) {
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
	public void printTSV(List<Feature> regionContents, OutputStream outputStream) {

		// Print header (row name column is nameless)
		PrintWriter out = new PrintWriter(outputStream);
		out.println("chromosome\tstart\tend\t" + Strings.delimit(regionContents.get(0).values.keySet(), "\t").toLowerCase());
		
		// Print row names and values
		HashSet<String> rowNames = new HashSet<String>();
		for (Feature regionContent : regionContents) {
			String rowName = regionContent.values.get(DataType.ID).toString();
			if (rowNames.contains(rowName)) {
				continue;
			}
			out.println(rowName + "\t" + regionContent.toString());
			rowNames.add(rowName);
		}
		
		out.flush();
	}

	public void sort(List<Feature> rows) {
		Collections.sort(rows);
	}
	
	/**
	 * Augment extra fields of primary RegionContent with secondary RegionContent.
	 * @return
	 */
	private static Feature augment(Feature primary, Feature secondary) {
		
		Feature augmented = new Feature(primary.region, getEmptyExtraFieldMap());
		augmented.values.clear();
		
		// Copy from primary, but augmenting from secondary
		for (Entry<DataType, Object> entry : primary.values.entrySet()) {
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
	

	public static LinkedHashMap<DataType, Object> getEmptyExtraFieldMap() {
		LinkedHashMap<DataType, Object> values = new LinkedHashMap<DataType, Object>();
		values.put(DataType.VALUE, EMPTY_EXTRA_FIELDS);
		return values;
	}
}
