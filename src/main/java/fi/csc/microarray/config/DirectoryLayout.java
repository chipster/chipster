package fi.csc.microarray.config;

import java.io.File;
import java.io.IOException;

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

	private static final String WORKSUBDIR = "nami-work-files";
	
	public enum Type {
		CLIENT,
		SERVER;
	}

	private Type type;
	private File workDir;
	private Configuration configuration;
	private static DirectoryLayout instance;

	public static DirectoryLayout initialiseServerLayout() throws IOException, OldConfigurationFormatException {
		synchronized (DirectoryLayout.class) {
			if (DirectoryLayout.instance != null) {
				throw new IllegalStateException("already initialised");
			}
			DirectoryLayout.instance = new DirectoryLayout(Type.SERVER, false, null);
			return DirectoryLayout.instance;
		}
	}

	public static DirectoryLayout initialiseClientLayout() throws IOException, OldConfigurationFormatException {
		return initialiseClientLayout(false, null);
	}

	public static DirectoryLayout initialiseClientLayout(boolean useHomeAsWorkDir, String overrideString) throws IOException, OldConfigurationFormatException {
		synchronized (DirectoryLayout.class) {
			if (DirectoryLayout.instance != null) {
				throw new IllegalStateException("already initialised");
			}
			DirectoryLayout.instance = new DirectoryLayout(Type.CLIENT, useHomeAsWorkDir, overrideString);
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
	
	private DirectoryLayout(Type type, boolean useHomeAsWorkDir, String overrideString) throws IOException, OldConfigurationFormatException {
		this.type = type;
		this.workDir = createWorkingDirectory(useHomeAsWorkDir);
		System.setProperty(Configuration.WORKDIR_PROPERTY[1], workDir.getAbsolutePath()); // NOTE: NO LOGGING IS TO BE DONE BEFORE THIS!
		this.configuration = new Configuration(overrideString, this.workDir.getAbsolutePath());
	}

	public File getWorkingDirectory() {
		return this.workDir;
	}
	
	private File createWorkingDirectory(boolean useHomeAsWorkDir) {
		String workDir;
		if (useHomeAsWorkDir) {
			String home = System.getProperty("user.home");
			workDir = home + File.separator + WORKSUBDIR;
		} else {
			workDir = new File(WORKSUBDIR).getAbsolutePath();
		}
		new File(workDir).mkdirs();

		return new File(workDir);
	}

	public File getFileroot() throws IOException, OldConfigurationFormatException {
		if (type == Type.SERVER) {
			File fileRepository = new File(getConfigurationValue("frontend", "fileServerPath"));
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
	
	private String getConfigurationValue(String module, String name) throws IOException, OldConfigurationFormatException {
		// we cannot assume config to be loaded inside this class (because configuration is managed via this class)
		getConfiguration();
		return Configuration.getValue(module, name);
	}

	public Configuration getConfiguration() throws IOException, OldConfigurationFormatException {
		return configuration;
	}	
}
