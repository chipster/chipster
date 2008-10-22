package fi.csc.microarray;

import java.io.IOException;
import java.util.Properties;

public class ApplicationConstants {
	
	static {
		try {
			Properties version = new Properties();
			version.load(ApplicationConstants.class.getResourceAsStream("/version.number"));
			Properties build = new Properties();
			build.load(ApplicationConstants.class.getResourceAsStream("/build.number"));
			NAMI_VERSION = version.getProperty("version.number") + " (build " + build.getProperty("build.number") + ")";

		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	public static final String NAMI_VERSION;
	public static final String APPLICATION_TITLE = "Chipster v" + NAMI_VERSION;
}
