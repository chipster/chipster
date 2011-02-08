package fi.csc.microarray.gbrowser.intervals;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.LinkedList;

import org.testng.Assert;
import org.testng.annotations.Test;

import fi.csc.microarray.client.gbrowser.intervals.IntervalOperations;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class IntervalToolTest {

	@Test
	public void test() throws FileNotFoundException, IOException {
		IntervalOperations tool = new IntervalOperations();

		LinkedList<RegionContent> rows1 = new LinkedList<RegionContent>();
		rows1.add(new RegionContent(new BpCoordRegion(100L, 200L, new Chromosome("1")), ""));
		rows1.add(new RegionContent(new BpCoordRegion(210L, 300L, new Chromosome("1")), ""));
		rows1.add(new RegionContent(new BpCoordRegion(400L, 500L, new Chromosome("1")), ""));
		rows1.add(new RegionContent(new BpCoordRegion(100L, 200L, new Chromosome("2")), ""));
		
		LinkedList<RegionContent> rows2 = new LinkedList<RegionContent>();
		rows2.add(new RegionContent(new BpCoordRegion(100L, 150L, new Chromosome("1")), ""));
		rows2.add(new RegionContent(new BpCoordRegion(250L, 600L, new Chromosome("1")), ""));
		rows2.add(new RegionContent(new BpCoordRegion(300L, 350L, new Chromosome("2")), ""));
		
		LinkedList<BpCoordRegion> expectedIntersection = new LinkedList<BpCoordRegion>();
		expectedIntersection.add(new BpCoordRegion(100L, 150L, new Chromosome("1")));
		expectedIntersection.add(new BpCoordRegion(250L, 300L, new Chromosome("1")));
		expectedIntersection.add(new BpCoordRegion(400L, 500L, new Chromosome("1")));
		Assert.assertEquals(tool.intersect(rows1, rows2, false), expectedIntersection);
		expectedIntersection.remove();
		Assert.assertNotSame(tool.merge(rows1, rows2), expectedIntersection);

		LinkedList<BpCoordRegion> expectedSubtraction = new LinkedList<BpCoordRegion>();
		expectedSubtraction.add(new BpCoordRegion(100L, 200L, new Chromosome("2")));
		Assert.assertEquals(tool.subtract(rows1, rows2), expectedSubtraction);
		expectedSubtraction.remove();
		Assert.assertNotSame(tool.merge(rows1, rows2), expectedSubtraction);

		LinkedList<BpCoordRegion> expectedMerge = new LinkedList<BpCoordRegion>();
		expectedMerge.add(new BpCoordRegion(100L, 200L, new Chromosome("1")));
		expectedMerge.add(new BpCoordRegion(210L, 600L, new Chromosome("1")));
		expectedMerge.add(new BpCoordRegion(100L, 200L, new Chromosome("2")));
		expectedMerge.add(new BpCoordRegion(300L, 350L, new Chromosome("2")));
		Assert.assertEquals(tool.merge(rows1, rows2), expectedMerge);
		expectedMerge.remove();
		Assert.assertNotSame(tool.merge(rows1, rows2), expectedMerge);

	}
	
	
	public static void main(String[] args) throws Exception {
		new IntervalToolTest().test();
		System.out.println("OK");
	}
}
