package fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat;

import org.testng.Assert;
import org.testng.annotations.Test;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.ChromosomeNameUnnormaliser;

public class ChromosomeNameUnnormaliserTest {

	@Test(groups = {"unit"} )
	public void test() {
		// Forward-reverse
		Assert.assertEquals(new ChromosomeNameUnnormaliser("X").unnormalise(new Chromosome("X")), "X");
		Assert.assertEquals(new ChromosomeNameUnnormaliser("chrX").unnormalise(new Chromosome("X")), "chrX");
		Assert.assertEquals(new ChromosomeNameUnnormaliser("chrX.fa").unnormalise(new Chromosome("X")), "chrX.fa");
		Assert.assertEquals(new ChromosomeNameUnnormaliser("X.fa").unnormalise(new Chromosome("X")), "X.fa");
		Assert.assertEquals(new ChromosomeNameUnnormaliser("chrX_random").unnormalise(new Chromosome("X_random")), "chrX_random");
		
		// Across naming conventions
		Assert.assertEquals(new ChromosomeNameUnnormaliser("X").unnormalise(new Chromosome("chrX")), "X");
		Assert.assertEquals(new ChromosomeNameUnnormaliser("chrX.fa").unnormalise(new Chromosome("X.fa")), "chrX.fa");
	}
	
	public static void main(String[] args) throws Exception {
		new ChromosomeNameUnnormaliserTest().test();
		System.out.println("PASSED");
	}
}
