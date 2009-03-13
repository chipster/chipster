package fi.csc.microarray.client;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;

public class ClientConstants {
	public final int MAX_JOBS;
	public final long MAX_JOB_SIZE;
	
	public ClientConstants() {
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();
		MAX_JOBS = Integer.parseInt(configuration.getValue("client", "max-jobs"));
		MAX_JOB_SIZE = Long.parseLong(configuration.getValue("client", "max-job-size"));
	}
}
