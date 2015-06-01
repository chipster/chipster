package fi.csc.microarray.jobmanager.model;

import java.util.Date;

public class Job {

	private int id; // primary key
	
	private String jobId;
	private String sessionId;
	private String description;
	private String headers;
	private String results;

	private String analysisId;
	private String compId;
	private String state;

	private String username;
	private String replyTo;
	
	private Date created;
	private Date rescheduled;
	private Date submitted;
	private Date finished;
	private Date seen;
	private Date dequeued;
	private Date explicitWait;
	
	private Integer retries;
	
	
}
