package fi.csc.microarray;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;

import javax.xml.parsers.ParserConfigurationException;

import org.mortbay.util.IO;
import org.xml.sax.SAXException;

import fi.csc.microarray.util.config.Configuration;
import fi.csc.microarray.util.config.ConfigurationLoader;
import fi.csc.microarray.util.config.ConfigurationLoader.OldConfigurationFormatException;

public class MicroarrayConfiguration {

//	private static final int REQUIRED_CONTENT_VERSION = 2;
	private static final String WORKSUBDIR = "nami-work-files";
	private static final File CONFIG_FILE = new File("nami-config.xml");	
	private static final String STATIC_CONFIG_RESOURCENAME = "/nami-static-config.xml";
	private static final String[] WORKDIR_PROPERTY = {null, "nami_work_dir"};	
	private static String DEFAULT_CONFIG_FILE = "/nami-config.xml.default";	
	private static boolean alreadyLoaded = false;	
	private static Configuration configuration = new Configuration(true);	
	
	/**
	 * Convenence method. Calls original wirh no broker override and without using home
	 * as work dir.
	 * @throws OldConfigurationFormatException 
	 * 
	 * @see #loadConfiguration(String, boolean)
	 */
	public static boolean loadConfiguration() throws IOException, OldConfigurationFormatException {
		return loadConfiguration(null, false);
	}
			
	/**
	 * <p>Initialises Nami-specific System-properties. If already loaded, then 
	 * does nothing, so calling this multiple times is safe. Method is thread
	 * safe.</p>
	 * 
	 * <p>If config file is not found, default settings are loaded and saved
	 * to config file for modifications.</p>
	 * 
	 * <p>NOTE: NO LOGGER IS TO BE CREATED BEFORE CALLING THIS METHOD!</p>
	 * 
	 * @throws IOException
	 * @return true if config was not found
	 * @throws OldConfigurationFormatException 
	 * @throws MicroarrayException 
	 */
	public static boolean loadConfiguration(String overrideString, boolean useHomeAsWorkDir) throws IOException, OldConfigurationFormatException {
		
		boolean defaultsUsed = false;
		
		// guard against reloading
		synchronized (MicroarrayConfiguration.class) {
			if (!alreadyLoaded) {

				// first set up work directory (for logging etc.)
				String workDir;
				if (useHomeAsWorkDir) {
					String home = System.getProperty("user.home");
					workDir = home + File.separator + WORKSUBDIR;
				} else {
					workDir = new File(WORKSUBDIR).getAbsolutePath();
				}
				new File(workDir).mkdirs();

				System.setProperty(WORKDIR_PROPERTY[1], workDir); // NOTE: NO LOGGING IS TO BE DONE BEFORE THIS!
				
				// if config file not available, copy defaults to it
				File configFile = new File(workDir + File.separator + CONFIG_FILE.getName());

				if (!configFile.exists()) {
					InputStream defaults = Configuration.class.getResourceAsStream(DEFAULT_CONFIG_FILE);
					OutputStream out = new FileOutputStream(configFile);
					
					try {
						IO.copy(defaults, out);
						defaultsUsed = true;
						
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
							configuration, MicroarrayConfiguration.class.getResourceAsStream(STATIC_CONFIG_RESOURCENAME), 0);
					ConfigurationLoader.addFromFile(configuration, configFile, 0);
					
				} catch (SAXException e) {
					throw new IOException(e.getMessage());
					
				} catch (ParserConfigurationException e) {
					throw new IOException(e.getMessage());
					
				} 

				// do overrides
				configuration.putValueInSubmodule(WORKDIR_PROPERTY[0], WORKDIR_PROPERTY[1], workDir);
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
			return defaultsUsed;
		}
	}

	public static String[] getValues(String name) {
		return getValues(null, name);
	}
	
	public static String[] getValues(String moduleName, String name) {
		assert(alreadyLoaded);
		Configuration module = getModule(moduleName);		
		return module.getValues(name);
	}

	public static String getValue(String name) {
		return getValue(null, name);
	}
	public static String getValue(String moduleName, String name) {
		assert(alreadyLoaded);
		Configuration module = getModule(moduleName); 
		return module.getValue(name);
	}

	private static Configuration getModule(String moduleName) {
		return moduleName != null ? configuration.getModule(moduleName) : configuration;
	}
	
	public static File getWorkDir() {
		return new File(getValue(WORKDIR_PROPERTY[0], WORKDIR_PROPERTY[1]));
	}
}
