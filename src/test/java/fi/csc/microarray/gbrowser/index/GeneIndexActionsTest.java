package fi.csc.microarray.gbrowser.index;

import java.sql.SQLException;

import org.testng.Assert;
import org.testng.annotations.Test;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;

public class GeneIndexActionsTest {

	private static class TestableGeneIndexActions extends GeneIndexActions {
		public TestableGeneIndexActions() throws ClassNotFoundException, SQLException {
			super(null, "test");
		}
	}

	@Test(groups = {"unit"} )
	public void test() throws SQLException, ClassNotFoundException {
		TestableGeneIndexActions tgia = new TestableGeneIndexActions();
		tgia.insertGene("1", "100", "200", "a2m");
		Assert.assertEquals(tgia.getLocation("a2m"), new BpCoordRegion(Long.parseLong("100"), Long.parseLong("200"), new Chromosome("1")));
	}
	
	public static void main(String[] args) throws SQLException, ClassNotFoundException {
		new GeneIndexActionsTest().test();
		System.out.println("ok");
	}
}
