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
	private Date rescheduled;
	private Date submitted;
	private Date finished;
	private Date seen;
	private Date dequeued;
	private Date explicitWait;

	//	private String analysisId;
	//	private String username;


	public Job(JobMessage jobMessage) {
		this.jobMessage = jobMessage;
		this.replyTo = jobMessage.getReplyTo();
		this.created = new Date();
		this.state = JobState.NEW;
	}

	public String getJobId() {
		return jobMessage.getJobId();
	}

	public JobMessage getJobMessage() {
		return jobMessage;
	}


	public Date getSubmitted() {
		return submitted;
	}

	public void setSubmitted(Date submitted) {
		this.submitted = submitted;
	}


	public ResultMessage getResults() {
		return results;
	}

	public void setResults(ResultMessage results) {
		this.results = results;
	}


	public JobState getState() {
		return state;
	}

	public void setState(JobState state) {
		this.state = state;
	}

	public Date getExplicitWait() {
		return explicitWait;
	}

	public void setExplicitWait(Date explicitWait) {
		this.explicitWait = explicitWait;
	}

	public String getCompId() {
		return compId;
	}

	public void setCompId(String compId) {
		this.compId = compId;
	}

	public Destination getReplyTo() {
		return replyTo;
	}

	public void setReplyTo(Destination replyTo) {
		this.replyTo = replyTo;
	}

	public long getSecondsSinceCreated() {
		return (System.currentTimeMillis() - created.getTime()) / 1000;
	}

	public long getSecondsSinceLastSeen() {
		return (System.currentTimeMillis() - seen.getTime()) / 1000;
	}

	public Date getCreated() {
		return created;
	}

	public void setCreated(Date created) {
		this.created = created;
	}

	public Date getRescheduled() {
		return rescheduled;
	}

	public void setRescheduled(Date rescheduled) {
		this.rescheduled = rescheduled;
	}

	public Date getFinished() {
		return finished;
	}

	public void setFinished(Date finished) {
		this.finished = finished;
	}

	public Date getSeen() {
		return seen;
	}

	public void setSeen(Date seen) {
		this.seen = seen;
	}

	public Date getDequeued() {
		return dequeued;
	}

	public void setDequeued(Date dequeued) {
		this.dequeued = dequeued;
	}

}
