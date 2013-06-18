package fi.csc.chipster.web.adminweb;

import java.io.IOException;
import java.util.Arrays;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;

public class ChipsterConfiguration {

	public static Configuration getConfiguration() throws IOException, IllegalConfigurationException {
		if (!DirectoryLayout.isInitialised()) {
			DirectoryLayout.initialiseServerLayout(Arrays.asList(new String[] {"manager"}));
		}

		return DirectoryLayout.getInstance().getConfiguration();
	}
	
	public static void init() throws IOException, IllegalConfigurationException {
		getConfiguration();
	}
}
