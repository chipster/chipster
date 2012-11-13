package fi.csc.microarray.manager.web;

import java.io.IOException;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;

public class ChipsterConfiguration {

	public static Configuration getConfiguration() {
		try {	
			if (DirectoryLayout.isInitialised()) {
				//already initialised by the Manager, run in same JVM
				//DirectoryLayout.initialiseServerLayout(Arrays.asList(new String[] {"manager"}));
			} else {
				
				// Not a real server, use any development server config (and show it's data)
				String configURL = "http://chipster-devel.csc.fi:8061/chipster-config.xml";
				//private final String configURL = "http://chipster.csc.fi/chipster-config.xml";
				//private final String configURL = "http://chipster.csc.fi/beta/chipster-config.xml";
				
				DirectoryLayout.initialiseSimpleLayout(configURL).getConfiguration();				
			}
		} catch (IOException e) {
			e.printStackTrace();
		} catch (IllegalConfigurationException e) {
			e.printStackTrace();
		}

    	return DirectoryLayout.getInstance().getConfiguration();
	}
	
	public static void init() {
		getConfiguration();
	}
}
