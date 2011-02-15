package fi.csc.microarray.gbrowser.intervals;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.LinkedList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import fi.csc.chipster.tools.gbrowser.regions.RegionOperations;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class RegionToolTest {

	@Test
	public void testSort() throws FileNotFoundException, IOException {
		RegionOperations tool = new RegionOperations();
		
		LinkedList<BpCoordRegion> rows = new LinkedList<BpCoordRegion>();
		rows.add(new BpCoordRegion(210L, 600L, new Chromosome("1")));
		rows.add(new BpCoordRegion(300L, 350L, new Chromosome("2")));
		rows.add(new BpCoordRegion(100L, 200L, new Chromosome("2")));
		rows.add(new BpCoordRegion(100L, 200L, new Chromosome("1")));

		LinkedList<BpCoordRegion> sortedRows = new LinkedList<BpCoordRegion>();
		sortedRows.add(new BpCoordRegion(100L, 200L, new Chromosome("1")));
		sortedRows.add(new BpCoordRegion(210L, 600L, new Chromosome("1")));
		sortedRows.add(new BpCoordRegion(100L, 200L, new Chromosome("2")));
		sortedRows.add(new BpCoordRegion(300L, 350L, new Chromosome("2")));

		tool.sort(rows);
		Assert.assertEquals(rows, sortedRows);
	}

	@Test
	public void test() throws FileNotFoundException, IOException {
		RegionOperations tool = new RegionOperations();

		// First (left) input is parsed and has track name & extra data (use chr prefix here)
		String fileContents = 
			"track name=pairedReads description=\"Clone Paired Reads\" useScore=1\n" +
			"chr1\t100\t200\tcloneA\t960\t+\t1000\t5000\t0\t2\t567,488\t0,3512\n" +
			"chr1\t210\t300\tcloneA\t960\t+\t1000\t5000\t0\t2\t567,488\t0,3512\n" +
			"chr1\t400\t500\tcloneA\t960\t+\t1000\t5000\t0\t2\t567,488\t0,3512\n" +
			"chr2\t100\t200\tcloneA\t960\t+\t1000\t5000\t0\t2\t567,488\t0,3512\n";
		List<RegionContent> rows1 = tool.parseString(fileContents);
		
		// Check that parsing was ok
		Assert.assertEquals(rows1.size(), 4);
		Assert.assertEquals(rows1.get(0).values.get(ColumnType.STRAND), Strand.FORWARD);
		Assert.assertEquals(rows1.get(0).values.get(ColumnType.ID), "cloneA");
		Assert.assertEquals(rows1.get(0).region.start.chr, new Chromosome("chr1.fa"));
		Assert.assertEquals(rows1.get(0).region.start.chr, new Chromosome("chr1"));
		Assert.assertEquals(rows1.get(0).region.start.chr, new Chromosome("1"));
		
		
		// Second (right) input given directly (don't use chr prefix here)
		LinkedList<RegionContent> rows2 = new LinkedList<RegionContent>();
		rows2.add(new RegionContent(new BpCoordRegion(100L, 150L, new Chromosome("1")), ""));
		rows2.add(new RegionContent(new BpCoordRegion(250L, 600L, new Chromosome("1")), ""));
		rows2.add(new RegionContent(new BpCoordRegion(300L, 350L, new Chromosome("2")), ""));
		
		// Test intersection
		LinkedList<BpCoordRegion> expectedIntersection = new LinkedList<BpCoordRegion>();
		expectedIntersection.add(new BpCoordRegion(100L, 150L, new Chromosome("1")));
		expectedIntersection.add(new BpCoordRegion(250L, 300L, new Chromosome("1")));
		expectedIntersection.add(new BpCoordRegion(400L, 500L, new Chromosome("1")));
		Assert.assertEquals(tool.intersect(rows1, rows2, 1L, RegionOperations.INTERSECT_PAIR_POLICY), expectedIntersection);
		expectedIntersection.remove();
		Assert.assertNotSame(tool.merge(rows1, rows2, 1L, true), expectedIntersection);

		// Test subtraction
		LinkedList<BpCoordRegion> expectedSubtraction = new LinkedList<BpCoordRegion>();
		expectedSubtraction.add(new BpCoordRegion(100L, 200L, new Chromosome("2")));
		Assert.assertEquals(tool.subtract(rows1, rows2, 1L), expectedSubtraction);
		expectedSubtraction.remove();
		Assert.assertNotSame(tool.merge(rows1, rows2, 1L, true), expectedSubtraction);

		// Test merge
		LinkedList<BpCoordRegion> expectedMerge = new LinkedList<BpCoordRegion>();
		expectedMerge.add(new BpCoordRegion(100L, 200L, new Chromosome("1")));
		expectedMerge.add(new BpCoordRegion(210L, 600L, new Chromosome("1")));
		expectedMerge.add(new BpCoordRegion(100L, 200L, new Chromosome("2")));
		expectedMerge.add(new BpCoordRegion(300L, 350L, new Chromosome("2")));
		Assert.assertEquals(tool.merge(rows1, rows2, 1L, true), expectedMerge);
		expectedMerge.remove();
		Assert.assertNotSame(tool.merge(rows1, rows2, 1L, true), expectedMerge);

		// Test min interlap size setting 
		LinkedList<BpCoordRegion> expectedSubtractionWithLongOverlap = new LinkedList<BpCoordRegion>();
		expectedSubtractionWithLongOverlap.add(new BpCoordRegion(100L, 200L, new Chromosome("1")));
		expectedSubtractionWithLongOverlap.add(new BpCoordRegion(210L, 300L, new Chromosome("1")));
		expectedSubtractionWithLongOverlap.add(new BpCoordRegion(100L, 200L, new Chromosome("2")));
		Assert.assertEquals(tool.subtract(rows1, rows2, 51L), expectedSubtractionWithLongOverlap);
		expectedSubtractionWithLongOverlap.remove();
		Assert.assertNotSame(tool.merge(rows1, rows2, 51L, true), expectedSubtractionWithLongOverlap);
	}
	
	
	public static void main(String[] args) throws Exception {
		new RegionToolTest().test();
		new RegionToolTest().testSort();
		System.out.println("OK");
	}
}
