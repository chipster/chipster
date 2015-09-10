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
	
	public Set<String> getWaitingJobs() {
		return waitingJobs;
	}


	public boolean updateJobScheduled(String jobId, String compId) {
		Job job = jobs.get(jobId);

		if (job == null) {
			logger.warn("update scheduled failed for non-existent job " + jobId);
			return false;
		}
		
		if (job.getFinished() != null) {
			logger.warn(String.format("cannot schedule an already finished job %s, state: %s", jobId, job.getState()));
			return false;
		}
		
		// update state
		waitingJobs.remove(job.getJobId());
		job.setState(JobState.SCHEDULED);
		
		job.setScheduled(new Date());
	    job.setCompId(compId);

	    return true;
	}
	
	
	/**
	 * 
	 * @param jobId
	 * @param state
	 * @param results null if no results
	 * @return
	 */
	public boolean updateJobFinished(String jobId, JobState state, ResultMessage results) {
		Job job = getJob(jobId);

		if (job == null) {
			logger.warn("update finished failed for non-existent job " + jobId);
			return false;
		}
		
		if (job.getFinished() != null) {
			logger.warn(String.format("cannot finish an already finished job %s, old state: %s, new state: %s", jobId, job.getState(), state));
			return false;
		}
		
		job.setFinished(new Date());
		job.setState(state);
		job.setResults(results);
		
		return true;
	}

	public boolean updateJobRunning(String jobId) {
		Job job = getJob(jobId);

		if (job == null) {
			logger.warn("update running failed for non-existent job " + jobId);
			return false;
		}
		
		if (job.getFinished() != null) {
			logger.warn("cannot put a finished job " + jobId + " to running state");
			return false;
		}

		job.setSeen(new Date());
		job.setState(JobState.RUNNING);
		return true;
	}
	
	public boolean updateJobWaiting(String jobId) {
		Job job = getJob(jobId);

		if (job == null) {
			logger.warn("update waiting failed for non-existent job " + jobId);
			return false;
		}
		
		if (job.getFinished() != null) {
			logger.warn(String.format("cannot put a finished job %s to wait", jobId));
			return false;
		}
		
		// TODO should this be denied?
		//if (job.getState() == JobState.SCHEDULED) {
		//	return false;
		//}
		
		job.setState(JobState.WAITING);
		waitingJobs.add(jobId);
		return true;
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


	/**
	 * 
	 * @param jobId
	 * @return true if the job can be cancelled
	 */
	public boolean updateJobCancelled(String jobId) {

		// remove from waiting if exists
		waitingJobs.remove(jobId);

		Job job = getJob(jobId);
		
		if (job == null) {
			return false;
		}
		
		// already finished
		if (job.getFinished() != null)  {
			return false;
		}
		
		// cancel
		job.setState(JobState.CANCELLED);
		job.setFinished(new Date());
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
