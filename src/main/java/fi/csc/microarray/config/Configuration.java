package fi.csc.microarray.config;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import javax.xml.parsers.ParserConfigurationException;

import org.mortbay.util.IO;
import org.xml.sax.SAXException;

import fi.csc.microarray.config.ConfigurationLoader.OldConfigurationFormatException;

public class Configuration {

	public static final String CONFIG_FILENAME = "chipster-config.xml";
	
	private static final String STATIC_CONFIG_RESOURCENAME = "/nami-static-config.xml";
	private static String DEFAULT_CONFIG_FILE = "/nami-config.xml.default";	
	private static boolean alreadyLoaded = false;	
	private static ConfigurationModule configuration = new ConfigurationModule(true);	
	
	public Configuration(String overrideString, File workDir) throws IOException, OldConfigurationFormatException {
		
		// guard against reloading
		synchronized (Configuration.class) {
			if (!alreadyLoaded) {

				// if config file not available, copy defaults to it
				File configFile = new File(workDir, CONFIG_FILENAME);

				if (!configFile.exists()) {
					InputStream defaults = ConfigurationModule.class.getResourceAsStream(DEFAULT_CONFIG_FILE);
					OutputStream out = new FileOutputStream(configFile);
					
					try {
						IO.copy(defaults, out);
						
					} finally {
						if (defaults != null) {
							defaults.close();
						}
						if (out != null) {
							out.close();
						}
					}			
				}

				// load configuration
				try {
					ConfigurationLoader.addFromStream(
							configuration, Configuration.class.getResourceAsStream(STATIC_CONFIG_RESOURCENAME), 0);
					ConfigurationLoader.addFromFile(configuration, configFile, 0);
					
				} catch (SAXException e) {
					throw new IOException(e.getMessage());
					
				} catch (ParserConfigurationException e) {
					throw new IOException(e.getMessage());
					
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
							configuration.putValuesInSubmodule(moduleName, entryName, values);
						} else {
							configuration.putValueInSubmodule(moduleName, entryName, values[0]);
						}
					}
					
				}
			}

			alreadyLoaded = true;
		}
	}

	public static String[] getValues(String name) {
		return getValues(null, name);
	}
	
	public static String[] getValues(String moduleName, String name) {
		assert(alreadyLoaded);
		ConfigurationModule module = getModule(moduleName);		
		return module.getValues(name);
	}

	public static String getValue(String name) {
		return getValue(null, name);
	}
	public static String getValue(String moduleName, String name) {
		assert(alreadyLoaded);
		ConfigurationModule module = getModule(moduleName); 
		return module.getValue(name);
	}

	private static ConfigurationModule getModule(String moduleName) {
		return moduleName != null ? configuration.getModule(moduleName) : configuration;
	}
}
