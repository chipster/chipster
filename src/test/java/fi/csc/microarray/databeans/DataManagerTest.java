package fi.csc.microarray.databeans;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import fi.csc.microarray.ModulesForTesting;
import fi.csc.microarray.config.MicroarrayConfiguration;
import fi.csc.microarray.config.ConfigurationLoader.OldConfigurationFormatException;
import fi.csc.microarray.databeans.fs.FSDataManager;

public class DataManagerTest {

	private DataManager manager; 
	
	@BeforeClass(alwaysRun = true)
	public void init() throws IOException, OldConfigurationFormatException {
		MicroarrayConfiguration.loadConfiguration();
		this.manager = new FSDataManager();
		ModulesForTesting.getModulesForTesting().plugFeatures(this.manager);
	}
	
	@Test(groups = {"unit"} )
	public void testDataTypes() throws IOException {
		File file = File.createTempFile("test", ".png");
		Assert.assertEquals(manager.guessContentType(file).getType(), "image/png");
	}
	
}
