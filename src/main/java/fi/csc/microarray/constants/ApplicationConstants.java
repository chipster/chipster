package fi.csc.microarray.constants;

import java.io.IOException;
import java.util.Properties;

public class ApplicationConstants {
	
	static {
		try {
			Properties version = new Properties();
			version.load(ApplicationConstants.class.getResourceAsStream("/version.number"));
			Properties build = new Properties();
			build.load(ApplicationConstants.class.getResourceAsStream("/build.number"));
			VERSION = version.getProperty("version.number") + " (build " + build.getProperty("build.number") + ")";

		} catch (IOException e) {
			throw new RuntimeException(e);
		}
	}
	public static final String VERSION;
	public static final String TITLE = "Chipster v" + VERSION;
}
