package fi.csc.microarray.databeans;

import java.io.File;
import java.io.IOException;

import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.Test;

import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.module.DefaultModules;

public class DataManagerTest {

	private DataManager manager; 
	
	@BeforeClass(alwaysRun = true)
	public void init() throws IOException, IllegalConfigurationException {
		DirectoryLayout.initialiseSimpleLayout().getConfiguration();			
		this.manager = new DataManager();
		DefaultModules.getDefaultModules().plugFeatures(this.manager);
	}
	
	@Test(groups = {"unit"} )
	public void testDataTypes() throws IOException {
		File file = File.createTempFile("test", ".png");
		Assert.assertEquals(manager.guessContentType(file).getType(), "image/png");
	}
	
}
