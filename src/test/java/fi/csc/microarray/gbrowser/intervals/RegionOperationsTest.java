package fi.csc.microarray.gbrowser.intervals;

import java.io.ByteArrayOutputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.util.LinkedHashMap;
import java.util.LinkedList;
import java.util.List;

import org.testng.Assert;
import org.testng.annotations.Test;

import fi.csc.chipster.tools.gbrowser.regions.RegionOperations;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class RegionOperationsTest {

	@Test
	public void testSort() throws FileNotFoundException, IOException {
		RegionOperations tool = new RegionOperations();
		
		LinkedList<RegionContent> rows = new LinkedList<RegionContent>();
		rows.add(new RegionContent(new Region(210L, 600L, new Chromosome("1")), ""));
		rows.add(new RegionContent(new Region(300L, 350L, new Chromosome("2")), ""));
		rows.add(new RegionContent(new Region(100L, 200L, new Chromosome("2")), ""));
		rows.add(new RegionContent(new Region(100L, 200L, new Chromosome("1")), ""));

		LinkedList<RegionContent> sortedRows = new LinkedList<RegionContent>();
		sortedRows.add(new RegionContent(new Region(100L, 200L, new Chromosome("1")), ""));
		sortedRows.add(new RegionContent(new Region(210L, 600L, new Chromosome("1")), ""));
		sortedRows.add(new RegionContent(new Region(100L, 200L, new Chromosome("2")), ""));
		sortedRows.add(new RegionContent(new Region(300L, 350L, new Chromosome("2")), ""));

		tool.sort(rows);
		Assert.assertEquals(rows, sortedRows);
	}

	@Test
	public void test() throws FileNotFoundException, IOException {
		RegionOperations tool = new RegionOperations();

		// First (left) input is parsed and has track name & extra data (use chr prefix here)
		String fileContents = 
			"chr1\t100\t200\tcloneA\t960\t+\t1000\t5000\t0\t2\t567,488\t0,3512\n" +
			"chr1\t210\t300\tcloneA\t960\t+\t1000\t5000\t0\t2\t567,488\t0,3512\n" +
			"chr1\t400\t500\tcloneA\t960\t+\t1000\t5000\t0\t2\t567,488\t0,3512\n" +
			"chr2\t100\t200\tcloneA\t960\t+\t1000\t5000\t0\t2\t567,488\t0,3512\n";
		String fileContentsWithHeader = 
			"track name=pairedReads description=\"Clone Paired Reads\" useScore=1\n" +
			fileContents;
		List<RegionContent> rows1 = tool.parseString(fileContentsWithHeader);
		
		// Check that parsing was ok
		Assert.assertEquals(rows1.size(), 4);
		Assert.assertEquals(rows1.get(0).values.get(ColumnType.STRAND), Strand.FORWARD);
		Assert.assertEquals(rows1.get(0).values.get(ColumnType.ID), "cloneA");
		Assert.assertEquals(rows1.get(0).region.start.chr, new Chromosome("chr1.fa"));
		Assert.assertEquals(rows1.get(0).region.start.chr, new Chromosome("chr1"));
		Assert.assertEquals(rows1.get(0).region.start.chr, new Chromosome("1"));
		
		// Check that track header row does not cause trouble
		Assert.assertEquals(tool.parseString(fileContents), rows1);
		
		// Check that printed output matches parsed input
		String normalisedFileContents = fileContents.replace("chr", "");
		OutputStream stringOut = new ByteArrayOutputStream();
		tool.print(rows1, stringOut);
		Assert.assertEquals(stringOut.toString(), normalisedFileContents);
		
		// Second (right) input given directly (don't use chr prefix here)
		LinkedList<RegionContent> rows2 = new LinkedList<RegionContent>();
		rows2.add(new RegionContent(new Region(100L, 150L, new Chromosome("1")), ""));
		rows2.add(new RegionContent(new Region(250L, 600L, new Chromosome("1")), ""));
		rows2.add(new RegionContent(new Region(300L, 350L, new Chromosome("2")), ""));
		
		// Test intersection
		LinkedList<RegionContent> expectedIntersection = new LinkedList<RegionContent>();
		expectedIntersection.add(new RegionContent(new Region(100L, 150L, new Chromosome("1")), ""));
		expectedIntersection.add(new RegionContent(new Region(250L, 300L, new Chromosome("1")), ""));
		expectedIntersection.add(new RegionContent(new Region(400L, 500L, new Chromosome("1")), ""));
		Assert.assertEquals(tool.intersect(rows1, rows2, 1L, RegionOperations.INTERSECT_PAIR_POLICY, true), expectedIntersection);
		expectedIntersection.remove();
		Assert.assertNotSame(tool.merge(rows1, rows2, 1L, true), expectedIntersection);

		// Test intersection while preserving originals
		LinkedList<RegionContent> expectedIntersection2 = new LinkedList<RegionContent>();
		LinkedHashMap<ColumnType, Object> values = new LinkedHashMap<ColumnType, Object>();
		values.put(ColumnType.ID, "cloneA");
		values.put(ColumnType.VALUE, 960L);
		values.put(ColumnType.STRAND, Strand.FORWARD);
		values.put(ColumnType.THICK_START, "1000");
		values.put(ColumnType.THICK_END, "5000");
		values.put(ColumnType.ITEM_RGB, "0");
		values.put(ColumnType.BLOCK_COUNT, "2");
		values.put(ColumnType.BLOCK_SIZES, "567,488");
		values.put(ColumnType.BLOCK_STARTS, "0,3512");
		expectedIntersection2.add(new RegionContent(new Region(100L, 200L, new Chromosome("1")), values));
		expectedIntersection2.add(new RegionContent(new Region(100L, 150L, new Chromosome("1")), ""));
		expectedIntersection2.add(new RegionContent(new Region(210L, 300L, new Chromosome("1")), values));
		expectedIntersection2.add(new RegionContent(new Region(250L, 600L, new Chromosome("1")), ""));
		expectedIntersection2.add(new RegionContent(new Region(400L, 500L, new Chromosome("1")), values));
		expectedIntersection2.add(new RegionContent(new Region(250L, 600L, new Chromosome("1")), ""));
		Assert.assertEquals(tool.intersect(rows1, rows2, 1L, RegionOperations.ORIGINALS_PAIR_POLICY, false), expectedIntersection2);
		expectedIntersection2.remove();
		Assert.assertNotSame(tool.merge(rows1, rows2, 1L, true), expectedIntersection2);

		// Test subtraction
		LinkedList<RegionContent> expectedSubtraction = new LinkedList<RegionContent>();
		expectedSubtraction.add(new RegionContent(new Region(100L, 200L, new Chromosome("2")), ""));
		Assert.assertEquals(tool.subtract(rows1, rows2, 1L), expectedSubtraction);
		expectedSubtraction.remove();
		Assert.assertNotSame(tool.merge(rows1, rows2, 1L, true), expectedSubtraction);

		// Test merge
		LinkedList<RegionContent> expectedMerge = new LinkedList<RegionContent>();
		expectedMerge.add(new RegionContent(new Region(100L, 200L, new Chromosome("1")), ""));
		expectedMerge.add(new RegionContent(new Region(210L, 600L, new Chromosome("1")), ""));
		expectedMerge.add(new RegionContent(new Region(100L, 200L, new Chromosome("2")), ""));
		expectedMerge.add(new RegionContent(new Region(300L, 350L, new Chromosome("2")), ""));
		Assert.assertEquals(tool.merge(rows1, rows2, 1L, true), expectedMerge);
		expectedMerge.remove();
		Assert.assertNotSame(tool.merge(rows1, rows2, 1L, true), expectedMerge);

		// Test min interlap size setting 
		LinkedList<RegionContent> expectedSubtractionWithLongOverlap = new LinkedList<RegionContent>();
		expectedSubtractionWithLongOverlap.add(new RegionContent(new Region(100L, 200L, new Chromosome("1")), ""));
		expectedSubtractionWithLongOverlap.add(new RegionContent(new Region(210L, 300L, new Chromosome("1")), ""));
		expectedSubtractionWithLongOverlap.add(new RegionContent(new Region(100L, 200L, new Chromosome("2")), ""));
		Assert.assertEquals(tool.subtract(rows1, rows2, 51L), expectedSubtractionWithLongOverlap);
		expectedSubtractionWithLongOverlap.remove();
		Assert.assertNotSame(tool.merge(rows1, rows2, 51L, true), expectedSubtractionWithLongOverlap);
	}
	
	
	public static void main(String[] args) throws Exception {
		new RegionOperationsTest().test();
		new RegionOperationsTest().testSort();
		System.out.println("OK");
	}
}
