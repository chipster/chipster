package fi.csc.microarray.jobmanager.model;

import java.util.Date;
import java.util.HashMap;
import java.util.LinkedHashSet;
import java.util.Map;
import java.util.Set;

import javax.jms.Destination;

import org.apache.log4j.Logger;

import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.messaging.message.ResultMessage;

public class JobManagerDB {

	private static Logger logger;
	private Map<String, Job> jobs;

	private LinkedHashSet<String> waitingJobs;

	public JobManagerDB() {
		logger = Logger.getLogger(JobManagerDB.class);
		jobs = new HashMap<String, Job>();
		waitingJobs = new LinkedHashSet<String>();
	}
	
	
	public void addJob(JobMessage jobMessage) {
		
		String jobId = jobMessage.getJobId();
		if (jobId == null || jobId.isEmpty()) {
			throw new NullPointerException("null or empty job id");
		}

		if (jobs.containsKey(jobId)) {
			throw new RuntimeException("job with id " + jobId + " already exists");
		}
		
		jobs.put(jobMessage.getJobId(), new Job(jobMessage));
		waitingJobs.add(jobId);
	}
	
	
	public Job getJob(String jobId) {
		return jobs.get(jobId);
	}
	
	public void updateJobSubmitted(String jobId, String compId) {
		Job job = jobs.get(jobId);

		// checks
		if (job == null) {
			throw new NullPointerException("job not found: " + jobId);
		}
		if (job.getResults() != null) {
			throw new IllegalStateException("job already finished: "+ jobId);
		}

		// update state
		waitingJobs.remove(job.getJobId());
		job.setState(JobState.SUBMITTED);
		
		job.setSubmitted(new Date());
	    job.setCompId(compId);
	}
	
	
	public void updateJobResults(String jobId, JobState state, ResultMessage results) {
		Job job = getJob(jobId);
		
		if (job.getCompId() == null) {
			logger.warn("adding results to job " + jobId + " with no comp id");
		}
		
		job.setFinished(new Date());
		job.setState(state);
		job.setResults(results);
	}

	public void updateJobRunning(String jobId) {
		Job job = getJob(jobId);
		
		if (job.getFinished() != null) {
			String s = "cannot put a finished job to running state";
			logger.error(s);
			throw new IllegalStateException(s);
		}

		job.setSeen(new Date());
		job.setState(JobState.RUNNING);
	}
	
	public void updateJobWaiting(String jobId) {
		Job job = getJob(jobId);
		
		if (job.getFinished() != null) {
			String s = "cannot put a finished job to wait";
			logger.error(s);
			throw new IllegalStateException(s);
		}
		
		if (job.getState() == JobState.SUBMITTED) {
			return;
		}
		
		job.setState(JobState.WAITING);
		job.setExplicitWait(new Date());
		waitingJobs.add(jobId);
	}


	/**
	 * 
	 * @param jobId
	 * @param newClientReplyTo
	 * @return null if job not found
	 */
	public Job updateJobReplyTo(String jobId, Destination newClientReplyTo) {
		Job job = getJob(jobId);
		
		if (job == null) {
			return null;
		}

		job.setReplyTo(newClientReplyTo);
		return job;
	}


	public boolean updateJobCancelled(String jobId) {

		// remove from waiting if exists
		waitingJobs.remove(jobId);

		Job job = getJob(jobId);
		
		if (job == null) {
			return false;
		}
		
		// already finished TODO check for state also?
		if (job.getFinished() != null)  {
			return false;
		}
		
		// cancel
		job.setFinished(new Date());
		job.setState(JobState.CANCELLED);
		return true;
	}

	public void updateJobError(String jobId) {
		// remove from waiting if exists
		waitingJobs.remove(jobId);

		Job job = getJob(jobId);
		
		if (job == null) {
			return;
		}
		
		job.setState(JobState.ERROR);;
	}

	
	
	public Set<String> getWaitingJobs() {
		return waitingJobs;
	}
//    job = session.query(Job).filter(Job.submitted == None).order_by(desc(Job.dequeued)).first()
//    if job:
//        job.dequeued = datetime.datetime.utcnow()
//    return job

	public void updateJobMaxWaitTimeReached(String jobId) {
		waitingJobs.remove(jobId);

		Job job = getJob(jobId);
		if (job == null) {
			return;
		}
		
		job.setState(JobState.EXPIRED_WAITING);
		job.setFinished(new Date());
	}
	
}
