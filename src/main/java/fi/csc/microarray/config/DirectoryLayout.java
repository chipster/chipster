package fi.csc.microarray.config;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Method;

import fi.csc.microarray.config.ConfigurationLoader.OldConfigurationFormatException;


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

	public static final String CONF_DIR = "conf";
	public static final String SECURITY_DIR = "security";
	public static final String LOGS_DIR = "logs";

	private static final String CONF_DIR_SYSTEM_PROPERTY = "chipster_conf_dir";
	private static final String LOGS_DIR_SYSTEM_PROPERTY = "chipster_logs_dir";
	private static final String SECURITY_DIR_SYSTEM_PROPERTY = "chipster_security_dir";

	public enum Type {
		CLIENT,
		SERVER;
	}

	private Type type;
	private Configuration configuration;
	private static DirectoryLayout instance;

	public static DirectoryLayout initialiseServerLayout() throws IOException, OldConfigurationFormatException {
		synchronized (DirectoryLayout.class) {
			if (DirectoryLayout.instance != null) {
				throw new IllegalStateException("already initialised");
			}
			DirectoryLayout.instance = new DirectoryLayout(Type.SERVER, null);
			return DirectoryLayout.instance;
		}
	}

	public static DirectoryLayout initialiseClientLayout() throws IOException, OldConfigurationFormatException {
		return initialiseClientLayout(null);
	}

	public static DirectoryLayout initialiseClientLayout(String overrideString) throws IOException, OldConfigurationFormatException {
		synchronized (DirectoryLayout.class) {
			if (DirectoryLayout.instance != null) {
				throw new IllegalStateException("already initialised");
			}
			DirectoryLayout.instance = new DirectoryLayout(Type.CLIENT, overrideString);
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
	
	private DirectoryLayout(Type type, String overrideString) throws IOException, OldConfigurationFormatException {
		this.type = type;
		System.setProperty(LOGS_DIR_SYSTEM_PROPERTY, getLogsDir().getAbsolutePath()); // NOTE: NO LOGGING IS TO BE DONE BEFORE THIS!
		System.setProperty(CONF_DIR_SYSTEM_PROPERTY, getConfDir().getAbsolutePath()); 
		System.setProperty(SECURITY_DIR_SYSTEM_PROPERTY, getSecurityDir().getAbsolutePath());
		this.configuration = new Configuration(overrideString, getConfDir());
	}

	public File getFileroot() throws IOException, OldConfigurationFormatException {
		if (type == Type.SERVER) {
			File fileRepository = new File(Configuration.getValue("frontend", "fileServerPath"));
			if (!fileRepository.exists()) {
				boolean ok = fileRepository.mkdir();
				if (!ok) {
					throw new IOException("could not create file root at " + fileRepository);
				}
			}
			return fileRepository;
			
		} else {
			throw new UnsupportedOperationException();
		}
	}
	
	public Configuration getConfiguration() throws IOException, OldConfigurationFormatException {
		return configuration;
	}	

	public File getConfDir() throws IOException {
		return initialise(new File(CONF_DIR));
	}

	public File getSecurityDir() throws IOException {
		return initialise(new File(SECURITY_DIR));
	}
	
	private File initialise(File dir) throws IOException {
		if (!dir.exists()) {
			boolean ok = dir.mkdirs(); // create whole path if does not exist 
			if (!ok) {
				throw new IOException("could not create directory path " + dir.getAbsolutePath());
			}
		}
		return dir;
	}

	private File getLogsDir() throws IOException {
		return initialise(new File(LOGS_DIR));
	}

	public File getUserDataDir() {
		return null;
	}
	
	public File getClientSettingsDir() throws IOException {
		String osName = System.getProperty("os.name");

		File dir = null;
		if (osName.startsWith("Mac OS")) {
			dir = null;

		} else if (osName.startsWith("Windows Vista")) {
			// %systemdrive%\ProgramData\Chipster
			dir = null;

		} else if (osName.startsWith("Windows")) {
			//%systemdrive%\Documents and Settings\All Users\Application Data\Chipster
			dir = null;

		} 
		if (dir != null) {
			try {
				dir.mkdirs();
			} catch (SecurityException se) {
				dir = null; // could not create, so can not use
			}
		}
		
		if (dir != null) {
			return dir;
			
		} else {
			// fall back to *nix kind of behaviour
			return initialise(new File(System.getProperty("user.home"), ".chipster"));
		}
	}


}
