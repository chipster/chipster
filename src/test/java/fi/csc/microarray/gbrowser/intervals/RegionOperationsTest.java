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

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.RegionOperations;

public class RegionOperationsTest {

	@Test
	public void testSort() throws FileNotFoundException, IOException {
		RegionOperations tool = new RegionOperations();
		
		LinkedList<RegionContent> rows = new LinkedList<RegionContent>();
		rows.add(new RegionContent(new Region(210L, 600L, new Chromosome("1")), RegionOperations.getEmptyExtraFieldMap()));
		rows.add(new RegionContent(new Region(300L, 350L, new Chromosome("2")), RegionOperations.getEmptyExtraFieldMap()));
		rows.add(new RegionContent(new Region(100L, 200L, new Chromosome("2")), RegionOperations.getEmptyExtraFieldMap()));
		rows.add(new RegionContent(new Region(100L, 200L, new Chromosome("1")), RegionOperations.getEmptyExtraFieldMap()));

		LinkedList<RegionContent> sortedRows = new LinkedList<RegionContent>();
		sortedRows.add(new RegionContent(new Region(100L, 200L, new Chromosome("1")), RegionOperations.getEmptyExtraFieldMap()));
		sortedRows.add(new RegionContent(new Region(210L, 600L, new Chromosome("1")), RegionOperations.getEmptyExtraFieldMap()));
		sortedRows.add(new RegionContent(new Region(100L, 200L, new Chromosome("2")), RegionOperations.getEmptyExtraFieldMap()));
		sortedRows.add(new RegionContent(new Region(300L, 350L, new Chromosome("2")), RegionOperations.getEmptyExtraFieldMap()));

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
		
		List<RegionContent> rows1 = new LinkedList<RegionContent>();
		for (String line : fileContentsWithHeader.split("\n")) {
			RegionContent region = tool.parseString(line);
			if (region != null) {
				rows1.add(region);
			}
		}
		
		// Check that parsing was ok
		Assert.assertEquals(rows1.size(), 4);
		Assert.assertEquals(rows1.get(0).values.get(DataType.STRAND), "+");
		Assert.assertEquals(rows1.get(0).values.get(DataType.ID), "cloneA");
		Assert.assertEquals(rows1.get(0).region.start.chr, new Chromosome("chr1.fa"));
		Assert.assertEquals(rows1.get(0).region.start.chr, new Chromosome("chr1"));
		Assert.assertEquals(rows1.get(0).region.start.chr, new Chromosome("1"));
		
		List<RegionContent> rowsNoHeader = new LinkedList<RegionContent>();
		for (String line : fileContents.split("\n")) {
			RegionContent region = tool.parseString(line);
			if (region != null) {
				rowsNoHeader.add(region);
			}
		}
		
		// Check that track header row does not cause trouble
		Assert.assertEquals(rowsNoHeader, rows1);
		
		// Check that printed output matches parsed input
		String normalisedFileContents = fileContents.replace("chr", "");
		OutputStream stringOut = new ByteArrayOutputStream();
		tool.print(rows1, stringOut);
		Assert.assertEquals(stringOut.toString(), normalisedFileContents);
		
		// Second (right) input given directly (don't use chr prefix here)
		LinkedList<RegionContent> rows2 = new LinkedList<RegionContent>();
		rows2.add(new RegionContent(new Region(100L, 150L, new Chromosome("1")), RegionOperations.getEmptyExtraFieldMap()));
		rows2.add(new RegionContent(new Region(250L, 600L, new Chromosome("1")), RegionOperations.getEmptyExtraFieldMap()));
		rows2.add(new RegionContent(new Region(300L, 350L, new Chromosome("2")), RegionOperations.getEmptyExtraFieldMap()));
		
		// Test intersection
		LinkedList<RegionContent> expectedIntersection = new LinkedList<RegionContent>();
		expectedIntersection.add(new RegionContent(new Region(100L, 150L, new Chromosome("1")), RegionOperations.getEmptyExtraFieldMap()));
		expectedIntersection.add(new RegionContent(new Region(250L, 300L, new Chromosome("1")), RegionOperations.getEmptyExtraFieldMap()));
		expectedIntersection.add(new RegionContent(new Region(400L, 500L, new Chromosome("1")), RegionOperations.getEmptyExtraFieldMap()));
		Assert.assertEquals(tool.intersect(rows1, rows2, 1L, RegionOperations.INTERSECT_PAIR_POLICY, true), expectedIntersection);
		expectedIntersection.remove();
		Assert.assertNotSame(tool.merge(rows1, rows2, 1L, true), expectedIntersection);

		// Test intersection while preserving originals
		LinkedList<RegionContent> expectedIntersection2 = new LinkedList<RegionContent>();
		LinkedHashMap<DataType, Object> values = new LinkedHashMap<DataType, Object>();
		values.put(DataType.ID, "cloneA");
		values.put(DataType.VALUE, 960L);
		values.put(DataType.STRAND, Strand.FORWARD);
		values.put(DataType.THICK_START, "1000");
		values.put(DataType.THICK_END, "5000");
		values.put(DataType.ITEM_RGB, "0");
		values.put(DataType.BLOCK_COUNT, "2");
		values.put(DataType.BLOCK_SIZES, "567,488");
		values.put(DataType.BLOCK_STARTS, "0,3512");
		expectedIntersection2.add(new RegionContent(new Region(100L, 200L, new Chromosome("1")), values));
		expectedIntersection2.add(new RegionContent(new Region(100L, 150L, new Chromosome("1")), RegionOperations.getEmptyExtraFieldMap()));
		expectedIntersection2.add(new RegionContent(new Region(210L, 300L, new Chromosome("1")), values));
		expectedIntersection2.add(new RegionContent(new Region(250L, 600L, new Chromosome("1")), RegionOperations.getEmptyExtraFieldMap()));
		expectedIntersection2.add(new RegionContent(new Region(400L, 500L, new Chromosome("1")), values));
		expectedIntersection2.add(new RegionContent(new Region(250L, 600L, new Chromosome("1")), RegionOperations.getEmptyExtraFieldMap()));
		Assert.assertEquals(tool.intersect(rows1, rows2, 1L, RegionOperations.ORIGINALS_PAIR_POLICY, false), expectedIntersection2);
		expectedIntersection2.remove();
		Assert.assertNotSame(tool.merge(rows1, rows2, 1L, true), expectedIntersection2);

		// Test subtraction
		LinkedList<RegionContent> expectedSubtraction = new LinkedList<RegionContent>();
		expectedSubtraction.add(new RegionContent(new Region(100L, 200L, new Chromosome("2")), RegionOperations.getEmptyExtraFieldMap()));
		Assert.assertEquals(tool.subtract(rows1, rows2, 1L), expectedSubtraction);
		expectedSubtraction.remove();
		Assert.assertNotSame(tool.merge(rows1, rows2, 1L, true), expectedSubtraction);

		// Test merge
		LinkedList<RegionContent> expectedMerge = new LinkedList<RegionContent>();
		expectedMerge.add(new RegionContent(new Region(100L, 200L, new Chromosome("1")), RegionOperations.getEmptyExtraFieldMap()));
		expectedMerge.add(new RegionContent(new Region(210L, 600L, new Chromosome("1")), RegionOperations.getEmptyExtraFieldMap()));
		expectedMerge.add(new RegionContent(new Region(100L, 200L, new Chromosome("2")), RegionOperations.getEmptyExtraFieldMap()));
		expectedMerge.add(new RegionContent(new Region(300L, 350L, new Chromosome("2")), RegionOperations.getEmptyExtraFieldMap()));
		Assert.assertEquals(tool.merge(rows1, rows2, 1L, true), expectedMerge);
		expectedMerge.remove();
		Assert.assertNotSame(tool.merge(rows1, rows2, 1L, true), expectedMerge);

		// Test min interlap size setting 
		LinkedList<RegionContent> expectedSubtractionWithLongOverlap = new LinkedList<RegionContent>();
		expectedSubtractionWithLongOverlap.add(new RegionContent(new Region(100L, 200L, new Chromosome("1")), RegionOperations.getEmptyExtraFieldMap()));
		expectedSubtractionWithLongOverlap.add(new RegionContent(new Region(210L, 300L, new Chromosome("1")), RegionOperations.getEmptyExtraFieldMap()));
		expectedSubtractionWithLongOverlap.add(new RegionContent(new Region(100L, 200L, new Chromosome("2")), RegionOperations.getEmptyExtraFieldMap()));
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
