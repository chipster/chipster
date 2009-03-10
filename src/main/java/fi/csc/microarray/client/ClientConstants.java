package fi.csc.microarray.client;

import fi.csc.microarray.config.Configuration;

public class ClientConstants {
	public static int MAX_JOBS = Integer.parseInt(Configuration.getValue("client", "max-jobs"));
	public static long MAX_JOB_SIZE = Long.parseLong(Configuration.getValue("client", "max-job-size"));
}
