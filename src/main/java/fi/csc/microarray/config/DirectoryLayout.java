package fi.csc.microarray.config;

import java.io.File;
import java.io.IOException;
import java.net.URL;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;


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
	public static final String CONF_DIR = "conf";

	public static final String WEB_ROOT = "web-root"; // TODO in future WEB_ROOT should be configurable (not easy because needs to be understood by Jetty)

	private static final String DEBUG_BASE_DIR = "debug-base-dir";

	private static final String CONF_DIR_SYSTEM_PROPERTY = "chipster_conf_dir";
	private static final String LOGS_DIR_SYSTEM_PROPERTY = "chipster_logs_dir";
	private static final String SECURITY_DIR_SYSTEM_PROPERTY = "chipster_security_dir";
	
	public enum Type {
		CLIENT,
		SERVER;
	}

	private Type type;
	private Configuration configuration = null;
	private boolean hasConfig;
	private static DirectoryLayout instance;

	public static DirectoryLayout initialiseServerLayout(List<String> specificModules) throws IOException, IllegalConfigurationException {
		synchronized (DirectoryLayout.class) {
			if (DirectoryLayout.instance != null) {
				throw new IllegalStateException("already initialised");
			}
			List<String> configModules = new LinkedList<String>(); 
			configModules.addAll(Arrays.asList(new String[] {"messaging", "security"}));
			configModules.addAll(specificModules);
			DirectoryLayout.instance = new DirectoryLayout(Type.SERVER, null, configModules, true);
			return DirectoryLayout.instance;
		}
	}

	public static DirectoryLayout initialiseUnitTestLayout() throws IOException, IllegalConfigurationException {
		synchronized (DirectoryLayout.class) {
			if (DirectoryLayout.instance != null) {
				throw new IllegalStateException("already initialised");
			}
			DirectoryLayout.instance = new DirectoryLayout(Type.CLIENT, null, Arrays.asList(new String[] {"messaging", "security"}), false);
			return DirectoryLayout.instance;
		}
	}

	public static DirectoryLayout initialiseClientLayout() throws IOException, IllegalConfigurationException {
		return initialiseClientLayout(null);
	}

	public static DirectoryLayout initialiseClientLayout(String configURL) throws IOException, IllegalConfigurationException {
		synchronized (DirectoryLayout.class) {
			if (DirectoryLayout.instance != null) {
				throw new IllegalStateException("already initialised");
			}
			DirectoryLayout.instance = new DirectoryLayout(Type.CLIENT, configURL, Arrays.asList(new String[] {"messaging", "security", "client"}), true);
			return DirectoryLayout.instance;
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
	
	private DirectoryLayout(Type type, String configURL, List<String> configModules, boolean hasConfig) throws IOException, IllegalConfigurationException {
		this.type = type;
		this.hasConfig = hasConfig;
		System.setProperty(LOGS_DIR_SYSTEM_PROPERTY, getLogsDir().getAbsolutePath()); // NOTE: NO LOGGING IS TO BE DONE BEFORE THIS!
		if (new File(getBaseDir(), SECURITY_DIR).exists()) {
			System.setProperty(SECURITY_DIR_SYSTEM_PROPERTY, getSecurityDir().getAbsolutePath());
		}
		if (hasConfig) {
			System.setProperty(CONF_DIR_SYSTEM_PROPERTY, getConfDir().getAbsolutePath()); 
			if (configURL == null) {
				this.configuration = new Configuration(getConfDir(), configModules);
			} else {
				this.configuration = new Configuration(new URL(configURL), configModules);
			}
		}
	}

	public File getConfDir() throws IOException {
		checkConfiguration();
		return check(new File(getBaseDir(), CONF_DIR));
	}


	public File getSecurityDir() throws IOException {
		return check(new File(getBaseDir(), SECURITY_DIR));
	}
	

	private File getLogsDir() throws IOException {
		return check(new File(getBaseDir(), LOGS_DIR));
	}

	public File getFileRoot() throws IOException, IllegalConfigurationException {
		if (type == Type.SERVER) {
			File fileRoot = new File(getBaseDir(), configuration.getString("filebroker", "file-root-path"));
			return initialise(fileRoot);
			
		} else {
			throw new UnsupportedOperationException();
		}
	}

	public File getCustomScriptsDir() throws IOException, IllegalConfigurationException {
		if (type == Type.SERVER) {
			File customScriptsDir = new File(getBaseDir(), configuration.getString("comp", "custom-scripts-dir"));
			return initialise(customScriptsDir);
			
		} else {
			throw new UnsupportedOperationException();
		}
	}

	public File getBackupDir() throws IOException, IllegalConfigurationException {
		if (type == Type.SERVER) {
			File backupDir = new File(getBaseDir(), configuration.getString("manager", "backup-dir"));
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
	
	private File getClientSettingsDir() throws IOException {
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
		if (type == Type.CLIENT) {
			// use OS specific dir for clients
			return getClientSettingsDir(); 
			
		} else {
			// use working dir as a base dir for server components
			File baseDir = new File(System.getProperty("user.dir"));
			
			// switch to debug dir if exists
			File debugRoot = new File(baseDir, DEBUG_BASE_DIR);
			if (debugRoot.exists()) {
				baseDir = debugRoot;
			}
			
			return baseDir;
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
		if (!hasConfig) {
			throw new IllegalStateException("directory layout has no configuration");
		}
	}
}
