package fi.csc.microarray.config;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;

public class Configuration {

	public static final String CONFIG_FILENAME = "chipster-config.xml";

	private static final int REQUIRED_CONFIGURATION_VERSION = 3;

	private static String CONFIG_SPECIFICATION_FILE = "/chipster-config-specification.xml";	
	private static ConfigurationModule rootModule = new ConfigurationModule();

	private List<String> configModules;	

	public Configuration(URL configUrl, List<String> configModules) throws IOException, IllegalConfigurationException {
		this(configUrl.openConnection().getInputStream(), configModules);
	}
	
	public Configuration(File workDir, List<String> configModules) throws IOException, IllegalConfigurationException {
		this(new FileInputStream(new File(workDir, CONFIG_FILENAME)), configModules);
	}
	
	public Configuration(InputStream configXml, List<String> configModules) throws IOException, IllegalConfigurationException {

		this.configModules = configModules;

		// load configuration specification and actual configuration XML		
		ConfigurationLoader loader = new ConfigurationLoader(this, REQUIRED_CONFIGURATION_VERSION);
		try {
			loader.addFromStream(Configuration.class.getResourceAsStream(CONFIG_SPECIFICATION_FILE), true);
			loader.addFromStream(configXml, false);

		} catch (SAXException e) {
			throw new IOException(e.getMessage());

		} catch (ParserConfigurationException e) {
			throw new IOException(e.getMessage());

		} 

		List<String> missingValues = rootModule.findMissingValues();
		if (!missingValues.isEmpty()) {
			throw new IllegalConfigurationException("configuration values missing: " + missingValues);
		}
	}

	public String[] getValues(String moduleName, String name) {
		ConfigurationModule module = getModule(moduleName);		
		return module.getValues(name);
	}

	public String getValue(String moduleName, String name) {
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
