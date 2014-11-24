package fi.csc.microarray.config;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.io.FilenameUtils;
import org.xml.sax.SAXException;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;

/**
 * <p>For accessing Chipster configuration.  The configuration is loaded in two steps. 
 * First an internal default configuration is loaded (chipster-config-specification.xml) 
 * and then the normal configuration file chipster-config.xml.</p>
 * 
 * <p>See <a href="http://chipster.wiki.sourceforge.net/ConfigurationSystem">http://chipster.wiki.sourceforge.net/ConfigurationSystem</a> for a more detailed description.</p>
 * 
 * @author Aleksi Kallio
 *
 */
public class Configuration {

	public static final String CONFIG_FILENAME = "chipster-config.xml";

	private static final int REQUIRED_CONFIGURATION_VERSION = 3;

	private static String CONFIG_SPECIFICATION_FILE = "/chipster-config-specification.xml";	
	private static ConfigurationModule rootModule = new ConfigurationModule();
	private URL configRootURL = null;

	private List<String> configModules;	

	public Configuration(URL configUrl, List<String> configModules) throws IOException, IllegalConfigurationException {
		this(configUrl.openConnection().getInputStream(), configModules);
		this.configRootURL = new URL(FilenameUtils.getFullPath(configUrl.toString()));
	}
	
	public Configuration(File workDir, List<String> configModules) throws IOException, IllegalConfigurationException {
		File configFile;
		String mostSpecificModuleName = configModules.size() > 0 ? configModules.get(configModules.size()-1) : null;

		if (new File(workDir, CONFIG_FILENAME).exists()) {
			configFile = new File(workDir, CONFIG_FILENAME); // normal config file
		} else if (new File(workDir, mostSpecificModuleName + "-" + CONFIG_FILENAME).exists()) {
			configFile = new File(workDir, mostSpecificModuleName + "-" + CONFIG_FILENAME); // component specific filename for debugging and development
		} else {
			throw new IOException("cannot find configuration file " + new File(workDir, CONFIG_FILENAME));
		}
		
		initialise(new FileInputStream(configFile), configModules);
		
		this.configRootURL = new File(FilenameUtils.getFullPath(configFile.getPath())).toURI().toURL();
	}

	public Configuration(List<String> configModules) throws IOException, IllegalConfigurationException {
		initialise((InputStream)null, configModules);
	}

	public Configuration(InputStream configXml, List<String> configModules) throws IOException, IllegalConfigurationException {
		initialise(configXml, configModules);
	}
	
	private void initialise(InputStream configXml, List<String> configModules) throws IOException, IllegalConfigurationException {
		
		this.configModules = configModules;

		// load configuration specification and actual configuration XML (if present)		
		ConfigurationLoader loader = new ConfigurationLoader(this, REQUIRED_CONFIGURATION_VERSION);
		try {
			loader.addFromStream(Configuration.class.getResourceAsStream(CONFIG_SPECIFICATION_FILE), true);
			if (configXml != null) {
				loader.addFromStream(configXml, false);
			}

		} catch (SAXException e) {
			throw new IOException(e.getMessage());

		} catch (ParserConfigurationException e) {
			throw new IOException(e.getMessage());

		} 

		if (configXml != null) {
			List<String> missingValues = rootModule.findMissingValues();
			if (!missingValues.isEmpty()) {
				throw new IllegalConfigurationException("configuration values missing: " + missingValues);
			}
		}
	}

	public String getString(String moduleName, String name) {
		ConfigurationModule module = findModule(moduleName); 
		return module.getEntry(name).getString();
	}

	public int getInt(String moduleName, String name) {
		ConfigurationModule module = findModule(moduleName); 
		return module.getEntry(name).getInt();
	}

	public boolean getBoolean(String moduleName, String name) {
		ConfigurationModule module = findModule(moduleName); 
		return module.getEntry(name).getBoolean();
	}
	
	private ConfigurationModule findModule(String moduleName) {
		return moduleName != null ? rootModule.getModule(moduleName) : rootModule;
	}
	
	public boolean isModuleEnabled(String moduleName) {
		return this.configModules.contains(moduleName);
	}

	public ConfigurationModule getRootModule() {
		return rootModule;
	}
	
	public URL getConfigRootURL() {
		return configRootURL;
	}
}
