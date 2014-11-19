package fi.csc.microarray.config;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import org.apache.log4j.PropertyConfigurator;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.filebroker.FileServer;


/**
 * <p>Specifies Chipster directory layout. This class if the decisive specification 
 * for directory layout. Directory layout is derived using three sources of 
 * information: 1) constant paths in this class, 2) detected runtime platform,
 * and 3) configuration. For a human readable description of directory layout see 
 * <a href="http://chipster.wiki.sourceforge.net/DirectoryLayout">http://chipster.wiki.sourceforge.net/DirectoryLayout</a>. 
 * </p>
 * 
 * <p>NOTE: NO LOGGER IS TO BE CREATED BEFORE INITIALISING DIRECTORYLAYOUT!</p>
 *  
 * @author Aleksi Kallio
 *
 */
public class DirectoryLayout {

	public static final String BIN_DIR = "bin";
	public static final String LOGS_DIR = "logs";
	public static final String SECURITY_DIR = "security";
	public static final String MODULES_DIR = "modules";
	public static final String CONF_DIR = "conf";
	public static final String LOCAL_ANNOTATION_DIR = "genomebrowser";
	
	public static final String WEB_ROOT = "web-root"; // TODO in future WEB_ROOT should be configurable (not easy because needs to be understood by Jetty)
	public static final String WEB_APPS_DIR = "webapps"; 

	private static final String CONF_DIR_SYSTEM_PROPERTY = "chipster_conf_dir";
	private static final String LOGS_DIR_SYSTEM_PROPERTY = "chipster_logs_dir";
	private static final String SECURITY_DIR_SYSTEM_PROPERTY = "chipster_security_dir";
	
	public enum Type {
		CLIENT,
		SERVER;
	}
	
	public enum AvailableConfiguration {
		NONE,
		DEFAULTS,
		FULL
	}

	private Type type;
	private Configuration configuration = null;
	private AvailableConfiguration availableConfiguration;
	private static File baseDirOverride = null;
	private static DirectoryLayout instance;
	
	/**
	 * Uninitialise (reset) directory layout. For testing purposes.
	 */
	public static void uninitialise() {
		instance = null;
	}
	
	/**
	 * Override default base dir with something else. For testing purposes.
	 */
	public static void setBaseDirOverride(File baseDirOverride) {
		DirectoryLayout.baseDirOverride = baseDirOverride;
	}
	
	public static DirectoryLayout initialiseServerLayout(List<String> specificModules)
	        throws IOException, IllegalConfigurationException {
	    return initialiseServerLayout(specificModules, null);
	}

	public static DirectoryLayout initialiseServerLayout(List<String> specificModules,
	        String configURL) throws IOException, IllegalConfigurationException {
		synchronized (DirectoryLayout.class) {
			if (DirectoryLayout.instance != null) {
				throw new IllegalStateException("already initialised");
			}
			List<String> configModules = new LinkedList<String>(); 
			configModules.addAll(Arrays.asList(new String[] {"messaging", "security"}));
			configModules.addAll(specificModules);
			DirectoryLayout.instance = new DirectoryLayout(Type.SERVER, configURL, configModules, AvailableConfiguration.FULL);
			return DirectoryLayout.instance;
		}
	}

	public static DirectoryLayout initialiseUnitTestLayout() throws IOException, IllegalConfigurationException {
		synchronized (DirectoryLayout.class) {
			if (DirectoryLayout.instance != null) {
				throw new IllegalStateException("already initialised");
			}
			DirectoryLayout.instance = new DirectoryLayout(Type.CLIENT, null, Arrays.asList(new String[] {"messaging", "security"}), AvailableConfiguration.NONE);
			return DirectoryLayout.instance;
		}
	}

	public static DirectoryLayout initialiseSimpleLayout() throws IOException, IllegalConfigurationException {
		return initialiseClientLayout(null);
	}

	public static DirectoryLayout initialiseSimpleLayout(String configURL) throws IOException, IllegalConfigurationException {
		return initialiseClientLayout(configURL);
	}

	public static DirectoryLayout initialiseClientLayout(String configURL) throws IOException, IllegalConfigurationException {
		synchronized (DirectoryLayout.class) {
			if (DirectoryLayout.instance != null) {
				throw new IllegalStateException("already initialised");
			}
			DirectoryLayout.instance = new DirectoryLayout(Type.CLIENT, configURL, Arrays.asList(new String[] {"messaging", "security", "client"}), AvailableConfiguration.FULL);
			return DirectoryLayout.instance;
		}
	}

	public static boolean isInitialised() {
		synchronized (DirectoryLayout.class) {
			return DirectoryLayout.instance != null;
		}
	}
		
	public static DirectoryLayout getInstance() {
		synchronized (DirectoryLayout.class) {
			if (DirectoryLayout.instance == null) {
				throw new IllegalStateException("not initialised");
			}
			return DirectoryLayout.instance;
		}
	}
	
	private DirectoryLayout(Type type, String configURL, List<String> configModules, AvailableConfiguration availableConfiguration) throws IOException, IllegalConfigurationException {
		this.type = type;
		this.availableConfiguration = availableConfiguration;
		
		System.setProperty(LOGS_DIR_SYSTEM_PROPERTY, getLogsDir().getAbsolutePath()); // NOTE: NO LOGGING IS TO BE DONE BEFORE THIS!
		
		// workaround for IcedTea-web (disable logging if OpenJDK)
		String javaRuntimeName = System.getProperty("java.runtime.name");
		if (type == Type.SERVER || javaRuntimeName == null || !javaRuntimeName.contains("OpenJDK")) {
			// enable logging for all server runtimes and all client runtimes that are not OpenJDK
			PropertyConfigurator.configure(getClass().getResourceAsStream("/log4j-enabled.properties")); // replaced with "enabled" config		
		}
		
		System.setProperty(SECURITY_DIR_SYSTEM_PROPERTY, getSecurityDir().getAbsolutePath());
		
		switch (availableConfiguration) {
			case DEFAULTS:
				this.configuration = new Configuration(configModules);
				break;
				
			case FULL:
				// check if config comes from file or from URL
				if (configURL == null) {
					System.setProperty(CONF_DIR_SYSTEM_PROPERTY, getConfDir().getAbsolutePath()); 
					this.configuration = new Configuration(getConfDir(), configModules);
				} else {
					this.configuration = new Configuration(new URL(configURL), configModules);
				}
				break;
				
			case NONE:
				// do nothing
				break;
		}
	}

	public File getConfDir() throws IOException {
		checkConfiguration();
		return check(new File(getBaseDir(), CONF_DIR));
	}

	public File getLocalAnnotationDir() throws IOException {
		return check(new File(getBaseDir(), LOCAL_ANNOTATION_DIR));
	}

	public File getSecurityDir() throws IOException {
		return check(new File(getBaseDir(), SECURITY_DIR));
	}	

	public File getModulesDir() throws IOException {
		return check(new File(getBaseDir(), MODULES_DIR));
	}	

	public File getLogsDir() throws IOException {
		return check(new File(getBaseDir(), LOGS_DIR));
	}

	public File getWebappsDir() throws IOException {
		if (type == Type.SERVER) {
			File webappsDir = new File(getBaseDir(), WEB_APPS_DIR);
			return initialise(webappsDir);
			
		} else {
			throw new UnsupportedOperationException();
		}
	}

	public File getFileRoot() throws IOException, IllegalConfigurationException {
		if (type == Type.SERVER) {
			File fileRoot = new File(getBaseDir(), configuration.getString("filebroker", "file-root-path"));
			File cacheRoot = new File(fileRoot, FileServer.CACHE_PATH);
			File storageRoot = new File(fileRoot, FileServer.STORAGE_PATH);
			File publicRoot = new File(fileRoot, configuration.getString("filebroker", "public-path"));
			initialise(fileRoot);
			initialise(cacheRoot);
			initialise(storageRoot);
			initialise(publicRoot);
			
			return fileRoot;
			
		} else {
			throw new UnsupportedOperationException();
		}
	}

	public File getManagerBackupDir() throws IOException, IllegalConfigurationException {
		if (type == Type.SERVER) {
			File backupDir = new File(getBaseDir(), configuration.getString("manager", "backup-dir"));
			return initialise(backupDir);
			
		} else {
			throw new UnsupportedOperationException();
		}
	}

	public File getFilebrokerMetadataBackupDir() throws IOException, IllegalConfigurationException {
		if (type == Type.SERVER) {
			File backupDir = new File(getBaseDir(), configuration.getString("filebroker", "metadata-backup-dir"));
			return initialise(backupDir);
			
		} else {
			throw new UnsupportedOperationException();
		}
	}

	
	public File getJobsDataDirBase(String id) throws IOException, IllegalConfigurationException {
		if (type == Type.SERVER) {
			File jobsDataDir = new File(getBaseDir(), configuration.getString("comp", "work-dir"));
			return initialise(new File(jobsDataDir, id));
			
		} else {
			throw new UnsupportedOperationException();
		}
	}

	public File getUserDataDir() throws IOException {
		if (type == Type.CLIENT) {

			File dir = null;
			File home = new File(System.getProperty("user.home"));

			if (new File(home, "My Documents").exists()) {
				dir = new File(home, "My Documents");

			} else if (new File(home, "Documents").exists()) {
				dir = new File(home, "Documents");

			} else {
				dir = home;
			}

			return dir;
			
		} else {
			throw new IllegalStateException("not supported for type SERVER");
		}
	}
	
	public File getClientSettingsDir() throws IOException {
		File dir = null;
		String osName = System.getProperty("os.name");

		// guess proper place for application settings
		if (osName.startsWith("Mac OS")) {
			dir = new File(System.getProperty("user.home"), "Library" + File.separator + "Application Support");

		} else if (osName.startsWith("Windows")) {
			dir = new File(System.getProperty("user.home"), "Local Settings" + File.separator + "Application Data");

		} 
		
		// initialise it
		if (dir != null) {
			dir = new File(dir, "Chipster");
			try {
				dir.mkdirs();
			} catch (SecurityException se) {
				dir = null; // could not create properly, so can not use this
			}
		}
		
		// if it did not work out, fall back to *nix mode 
		if (dir == null) {
			dir = new File(System.getProperty("user.home"), ".chipster");
		}

		return check(dir);
	}

	private File getBaseDir() throws IOException {

		if (baseDirOverride != null) {
			return baseDirOverride;
			
		} else if (type == Type.CLIENT) {
			// use OS specific dir for clients
			return getClientSettingsDir(); 
			
		} else {
			// use working dir as a base dir for server components
			return new File(System.getProperty("user.dir"));
		}
	}
	
	private File check(File dir) throws IOException {
		if (!dir.exists()) {
			if (type == Type.CLIENT) {
				initialise(dir);
				
			} else {
				throw new IOException("directory " + dir.getAbsolutePath() + " does not exist");
			}
		}
		return dir;
	}
	
	private File initialise(File dir) throws IOException {
		dir.mkdirs(); // create whole path if does not exist 
		return dir;
	}

	public Configuration getConfiguration() {
		checkConfiguration();
		return configuration;
	}	

	private void checkConfiguration() {
		if (availableConfiguration == AvailableConfiguration.NONE) {
			throw new IllegalStateException("directory layout has no configuration");
		}
	}
	
	public Type getType() {
		return type;
	}
}
