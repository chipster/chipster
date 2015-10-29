package fi.csc.microarray.jobmanager.model;

import java.util.ArrayList;
import java.util.Date;
import java.util.List;
import java.util.Timer;
import java.util.TimerTask;
import java.util.UUID;

import javax.jms.Destination;

import org.apache.activemq.command.ActiveMQTempTopic;
import org.apache.log4j.Logger;
import org.joda.time.DateTime;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.jobmanager.HibernateUtil;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.messaging.message.ResultMessage;

/**
 * 
 * To access this database on the command line with default settings:
 * 
 * cd ~/workspace/chipster-environment
 * java -cp ../chipster/ext/lib/h2-1.3.163.jar org.h2.tools.Shell -url jdbc:h2:database/jobmanager-db -user sa -password ""
 * 
 * @author klemela
 *
 */
public class JobManagerDB {	

	private static Logger logger;
	private HibernateUtil hibernate;
	private Timer purgeOldJobsTimer;
	
	private int purgeJobsOlderThan;	
	private long purgeOldJobsInterval = 0;  
	
	public class PurgeOldJobsTask extends TimerTask {
		@Override
		public void run() {
			purgeOldJobs();
		}	
	}

	public JobManagerDB(Configuration configuration) {
		logger = Logger.getLogger(JobManagerDB.class);

		List<Class<?>> hibernateClasses = new ArrayList<Class<?>>();
		hibernateClasses.add(Job.class);
		this.hibernate = new HibernateUtil();
		this.hibernate.buildSessionFactory(hibernateClasses, configuration);
		
		this.purgeOldJobsInterval = configuration.getInt("jobmanager", "purge-jobs-interval"); // hours
		this.purgeJobsOlderThan = configuration.getInt("jobmanager", "purge-jobs-older-than"); // days
		
		logger.info("check for old jobs every " + purgeOldJobsInterval + " hours and purge all jobs older than " + purgeJobsOlderThan + " days");
		logger.info("there are " + getJobCount() + " jobs in the database");
		
		if (purgeOldJobsInterval > 0) {
			purgeOldJobsTimer = new Timer(true);
			purgeOldJobsTimer.schedule(new PurgeOldJobsTask(), purgeOldJobsInterval * 60l * 60 * 1000, purgeOldJobsInterval * 60l * 60 * 1000);
			
			logger.info("check for old jobs every " + purgeOldJobsInterval + " hours and purge all jobs older than " + purgeJobsOlderThan + " days");			
		} else {
			logger.info("check for old jobs is disabled");
		}
		logger.info("there are " + getJobCount() + " jobs in the database");
	}


	public boolean addJob(JobMessage jobMessage) {

		String jobId = jobMessage.getJobId();
		if (jobId == null || jobId.isEmpty()) {
			logger.warn("add job failed, jobId is: " + jobId);
			return false;
		}

		if (getJob(jobId) != null) {
			logger.warn("add job failed, job with id " + jobId + " already exists");
			return false;
		}

		Job job = new Job(jobMessage);
		job.setReplyTo(jobMessage.getReplyTo());

		this.hibernate.beginTransaction();
		try {
			this.hibernate.session().save(job);
			this.hibernate.commit();

			return true;
		} catch (Throwable e) {
			this.hibernate.rollback();
			throw e;
		}
	}


	public Job getJob(String jobId) {
		this.hibernate.beginTransaction();
		try {
			Job job = (Job) this.hibernate.session().get(Job.class, UUID.fromString(jobId));
			this.hibernate.commit();
			return job;
		} catch (Throwable e) {
			this.hibernate.rollback();
			throw e;
		}
	}

	public Job updateJob(Job job) {
		this.hibernate.beginTransaction();
		try {
			this.hibernate.session().merge(job);
			this.hibernate.commit();
			return job;
		} catch (Throwable e) {
			this.hibernate.rollback();
			throw e;
		}
	}

	public List<Job> getWaitingJobs() {
		this.hibernate.beginTransaction();
		try {
			@SuppressWarnings("unchecked")
			List<Job> jobs = this.hibernate.session().createQuery(
					"from Job "
							+ "where state=:state "
							+ "order by created")
							.setParameter("state", JobState.WAITING).list();

			this.hibernate.commit();
			return jobs;
		} catch (Throwable e) {
			this.hibernate.rollback();
			throw e;
		}
	}

	public List<Job> getRunningJobs() {
		this.hibernate.beginTransaction();
		try {
			@SuppressWarnings("unchecked")
			List<Job> jobs = this.hibernate.session().createQuery(
					"from Job "
							+ "where state=:state1 "
							+ "or state=:state2 "
							+ "order by created")
							.setParameter("state1", JobState.WAITING)
							.setParameter("state2", JobState.RUNNING)
							.list();			

			this.hibernate.commit();
			return jobs;
		} catch (Throwable e) {
			this.hibernate.rollback();
			throw e;
		}
	}


	public boolean updateJobScheduled(Job job, String compId, String compHost) {
		if (job == null) {
			logger.warn("update scheduled failed: job is null");
			return false;
		}

		if (job.getFinished() != null) {
			logger.warn(String.format("cannot schedule an already finished job %s, state: %s", job.getJobId(), job.getState()));
			return false;
		}

		// update state
		job.setState(JobState.SCHEDULED);

		job.setScheduled(new Date());
		job.setCompId(compId);
		job.setCompHost(compHost);

		updateJob(job);

		return true;
	}


	/**
	 * 
	 * @param job
	 * @param state
	 * @param results null if no results
	 * @return
	 */
	public boolean updateJobFinished(Job job, JobState state, ResultMessage results) {
		if (job == null) {
			logger.warn("update finished failed: job is null");
			return false;
		}

		if (job.getFinished() != null) {
			logger.warn(String.format("cannot finish an already finished job %s, old state: %s, new state: %s", job.getJobId(), job.getState(), state));
			return false;
		}

		job.setFinished(new Date());
		job.setState(state);
		job.setResults(results);

		updateJob(job);

		return true;
	}

	public boolean updateJobRunning(Job job) {
		if (job == null) {
			logger.warn("update running failed: job is null");
			return false;
		}

		if (job.getFinished() != null) {
			logger.warn("cannot put a finished job " + job.getJobId() + " to running state");
			return false;
		}

		job.setSeen(new Date());
		job.setState(JobState.RUNNING);

		updateJob(job);

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

		job.setReplyTo((ActiveMQTempTopic) newClientReplyTo);

		updateJob(job);

		return job;
	}


	/**
	 * 
	 * @param jobId
	 * @return true if the job can be cancelled
	 */
	public boolean updateJobCancelled(Job job) {
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

		updateJob(job);

		return true;
	}

	public void updateJobMaxWaitTimeReached(String jobId) {

		Job job = getJob(jobId);
		if (job == null) {
			return;
		}

		job.setState(JobState.EXPIRED_WAITING);
		job.setFinished(new Date());

		updateJob(job);
	}


	/**
	 * Not used at the moment
	 * 
	 * @param jobId
	 */
	public void updateJobError(String jobId) {

		Job job = getJob(jobId);

		if (job == null) {
			return;
		}

		job.setState(JobState.ERROR);;

		updateJob(job);
	}


	/**
	 * Not used at the moment
	 * 
	 * @param jobId
	 * @return
	 */
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

		updateJob(job);
		return true;
	}
	
	public void purgeOldJobs() {
		this.hibernate.beginTransaction();
		try {
			
			Date purgeDate = new DateTime().minusDays(purgeJobsOlderThan).toDate();

			logger.info("removing jobs older than " + purgeDate);
			int count = this.hibernate.session().createQuery(
					"delete from Job "
					+ "where created < :date ")
					.setParameter("date", purgeDate).executeUpdate();
			
			this.hibernate.commit();
			logger.info("jobs removed: " + count);

		} catch (Throwable e) {
			this.hibernate.rollback();
			throw e;
		}
	}
	
	public Long getJobCount() {
		this.hibernate.beginTransaction();
		try {			
			Long count = (Long) this.hibernate.session().createQuery(
					"select count(*) from Job").uniqueResult();			
			this.hibernate.commit();
			return count;
			
		} catch (Throwable e) {
			this.hibernate.rollback();
			throw e;
		}
	}	
}
