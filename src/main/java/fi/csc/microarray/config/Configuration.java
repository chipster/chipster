package fi.csc.microarray.config;

import java.io.File;
import java.io.IOException;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;

public class Configuration {

	public static final String CONFIG_FILENAME = "chipster-config.xml";

	private static final int REQUIRED_CONFIGURATION_VERSION = 3;
	
	private static String CONFIG_SPECIFICATION_FILE = "/chipster-config-specification.xml";	
	private static boolean alreadyLoaded = false;	
	private static ConfigurationModule rootModule = new ConfigurationModule();

	private List<String> configModules;	
	
	public Configuration(String overrideString, File workDir, List<String> configModules) throws IOException, IllegalConfigurationException {
		
		// guard against reloading
		synchronized (Configuration.class) {
			if (!alreadyLoaded) {
				
				this.configModules = configModules;

				// load configuration spec. 
				File configFile = new File(workDir, CONFIG_FILENAME);
				ConfigurationLoader loader = new ConfigurationLoader(this, REQUIRED_CONFIGURATION_VERSION);
				try {
					loader.addFromStream(Configuration.class.getResourceAsStream(CONFIG_SPECIFICATION_FILE), true);
					loader.addFromFile(configFile, false);
					
				} catch (SAXException e) {
					throw new IOException(e.getMessage());
					
				} catch (ParserConfigurationException e) {
					throw new IOException(e.getMessage());
					
				} 

				List<String> missingValues = rootModule.findMissingValues();
				if (!missingValues.isEmpty()) {
					throw new IllegalConfigurationException("configuration values missing: " + missingValues);
				}
				
				// do overrides
				if (overrideString != null) {
					String[] overrides = overrideString.split(",");
					for (String override : overrides) {
						// split to name and value parts
						String[] parts = override.split("=");
						
						// parse entry name
						String[] nameParts = parts[0].split("/");
						String moduleName = nameParts.length > 1 ? nameParts[0] : null;
						String entryName = nameParts.length > 1 ? nameParts[1] : nameParts[0];
						
						// parse entry value(s)
						String values[] = parts[1].split(";");
						if (values.length > 1) {
							rootModule.putValuesInSubmodule(moduleName, entryName, values);
						} else {
							rootModule.putValueInSubmodule(moduleName, entryName, values[0]);
						}
					}
					
				}
			}

			alreadyLoaded = true;
		}
	}

	public static String[] getValues(String moduleName, String name) {
		assert(alreadyLoaded);
		ConfigurationModule module = getModule(moduleName);		
		return module.getValues(name);
	}

	public static String getValue(String moduleName, String name) {
		assert(alreadyLoaded);
		ConfigurationModule module = getModule(moduleName); 
		return module.getValue(name);
	}

	private static ConfigurationModule getModule(String moduleName) {
		return moduleName != null ? rootModule.getModule(moduleName) : rootModule;
	}
	
	public boolean isModuleEnabled(String moduleName) {
		return this.configModules.contains(moduleName);
	}

	public ConfigurationModule getRootModule() {
		return rootModule;
	}

}
