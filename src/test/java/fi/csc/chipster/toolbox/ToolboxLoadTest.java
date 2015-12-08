package fi.csc.chipster.toolbox;

import java.io.File;
import java.io.IOException;

import org.junit.Test;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.config.DirectoryLayout;

public class ToolboxLoadTest {

	
	@Test
	public void loadToolbox() throws IOException, IllegalConfigurationException {
		DirectoryLayout.uninitialise();
		DirectoryLayout.initialiseUnitTestLayout();			
		
		@SuppressWarnings("unused")
		Toolbox toolbox = new Toolbox(new File("../chipster-tools/modules"));
	}
}
