package fi.csc.microarray.client;

import fi.csc.microarray.MicroarrayConfiguration;

public class ClientConstants {
	public static int MAX_JOBS = Integer.parseInt(MicroarrayConfiguration.getValue("client", "max_jobs"));
	public static long MAX_JOB_SIZE = Long.parseLong(MicroarrayConfiguration.getValue("client", "max_job_size"));
}
