package fi.csc.microarray.client;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;

public class ClientConstants {
	public final int MAX_JOBS;
	
	public ClientConstants() {
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
		MAX_JOBS = configuration.getInt("client", "max-jobs");
	}
}
