package fi.csc.microarray.client.gbrowser.intervals;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Arrays;
import java.util.Collections;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Chunk;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.BEDParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class IntervalTool {

	
	public void initialise() throws IOException {
		File input1 = new File("src/test/java/fi/csc/microarray/gbrowser/intervals/input1.txt");
		File input2 = new File("src/test/java/fi/csc/microarray/gbrowser/intervals/input2.txt");
		
		List<RegionContent> rows1 = loadFile(input1);
		List<RegionContent> rows2 = loadFile(input2);

		LinkedList<BpCoordRegion> intersections = intersect(rows1, rows2);
		print(intersections);
	}

	/**
	 * For each interval in first set, finds intersection with second set.
	 * Returns merged intersections.
	 *  
	 * @param rows1 first set 
	 * @param rows2 second set
	 * @return intersections
	 */
	private LinkedList<BpCoordRegion> intersect(List<RegionContent> rows1, List<RegionContent> rows2) {
		LinkedList<BpCoordRegion> intersections = new LinkedList<BpCoordRegion>();
		for (RegionContent row1 : rows1) {
			for (RegionContent row2 : rows2) {
				if (row1.region.intercepts(row2.region)) {
					intersections.add(row1.region.intercept(row2.region));
				}
			}
		}
		LinkedList<BpCoordRegion> mergedIntersections = mergeContinuous(intersections);
		return mergedIntersections;
	}

	private void print(LinkedList<BpCoordRegion> intervals) {
		for (BpCoordRegion interval : intervals) {
			System.out.println(interval);
		}
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
				if (intervals.get(i).intercepts(intervals.get(j + 1))) {
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

	private List<RegionContent> loadFile(File input1) throws FileNotFoundException, IOException {
		ChunkDataSource dataSource = new ChunkDataSource(input1, new BEDParser());
		byte[] fileChunk = dataSource.readAll();
		List<ColumnType> columns = Arrays.asList(new ColumnType[] { ColumnType.CHROMOSOME, ColumnType.BP_START, ColumnType.BP_END});
		return dataSource.getFileParser().getAll(new Chunk(new String(fileChunk)), columns);
	}
	
	public static void main(String[] args) throws IOException {
		new IntervalTool().initialise();
	}
}
