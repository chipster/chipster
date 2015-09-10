package fi.csc.microarray.jobmanager.model;

import java.util.Date;

import javax.jms.Destination;

import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.messaging.message.ResultMessage;

public class Job {

	private JobMessage jobMessage;
	private ResultMessage results;
	private String compId;
	private JobState state;
	private Destination replyTo;

	private Date created;
	private Date scheduled;
	private Date finished;
	private Date seen;

//	private Date rescheduled;
//	private Date dequeued;
//	private Date explicitWait;
//	private String analysisId;
//	private String username;


	Job(JobMessage jobMessage) {
		this.jobMessage = jobMessage;
		this.replyTo = jobMessage.getReplyTo();
		this.created = new Date();
		this.state = JobState.WAITING;
	}

	public String getJobId() {
		return jobMessage.getJobId();
	}

	public JobMessage getJobMessage() {
		return jobMessage;
	}

	public Date getScheduled() {
		return scheduled;
	}

	public ResultMessage getResults() {
		return results;
	}

	public JobState getState() {
		return state;
	}

	public String getCompId() {
		return compId;
	}

	public Destination getReplyTo() {
		return replyTo;
	}

	public Date getCreated() {
		return created;
	}

	public Date getFinished() {
		return finished;
	}

	public void setFinished(Date finished) {
		this.finished = finished;
	}

	public long getSecondsSinceCreated() {
		return (System.currentTimeMillis() - created.getTime()) / 1000;
	}

	void setScheduled(Date scheduled) {
		this.scheduled = scheduled;
	}

	void setResults(ResultMessage results) {
		this.results = results;
	}

	void setState(JobState state) {
		this.state = state;
	}

	void setCompId(String compId) {
		this.compId = compId;
	}

	void setReplyTo(Destination replyTo) {
		this.replyTo = replyTo;
	}

	void setSeen(Date seen) {
		this.seen = seen;
	}

}
