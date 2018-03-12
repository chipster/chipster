package fi.csc.microarray.jobmanager;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedList;
import java.util.List;

import javax.jms.Destination;
import javax.jms.JMSException;

import org.apache.log4j.Logger;

import com.google.gson.Gson;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.jobmanager.model.Job;
import fi.csc.microarray.jobmanager.model.JobManagerDB;
import fi.csc.microarray.messaging.JMSMessagingEndpoint;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingListener;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.MonitoredNodeBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.JobLogMessage;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.messaging.message.JsonMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;
import fi.csc.microarray.messaging.message.ResultMessage;
import fi.csc.microarray.messaging.message.ServerStatusMessage;
import fi.csc.microarray.service.KeepAliveShutdownHandler;
import fi.csc.microarray.service.ShutdownCallback;
import fi.csc.microarray.util.Exceptions;
import fi.csc.microarray.util.SystemMonitorUtil;

public class JobManager extends MonitoredNodeBase implements MessagingListener, ShutdownCallback {

	private int jobMaxWaitTime;
	
	private static Logger logger;

	private MessagingEndpoint endpoint;
	
	MessagingTopic fromClientTopic;
	MessagingTopic jobManagerAdminTopic;	
	MessagingTopic jobManagerTopic;
	MessagingTopic compTopic;
	
	private JobManagerDB jobsDb;
	
	
	private class ClientMessageListener implements MessagingListener {

		@Override
		public void onChipsterMessage(ChipsterMessage msg) {			
			try {
				
				// JobMessage
				if (msg instanceof JobMessage) {
					handleClientJobMessage((JobMessage) msg);
				} 

				// CommandMessage
				else if (msg instanceof CommandMessage) {
					handleClientCommandMessage((CommandMessage) msg);
				}

				// other messages
				else {	
					logger.warn("sending " + msg.getClass() + " directly to comps, replyTo is unchanged");
					// msg.setReplyTo(jobManagerTopic.getJMSTopic());
					compTopic.sendMessage(msg);
				}
			} catch (Exception e) {
				logger.error(Exceptions.getStackTrace(e));
			}
		}
	
	
		private void handleClientJobMessage(JobMessage msg) throws JMSException {
			
			Destination clientReplyTo = msg.getReplyTo();

			// store in db, also stores the original replyTo
			// when created, job is put to state waiting
			if (jobsDb.addJob(msg)) {  

				// set replyTo to jobmanager
				msg.setReplyTo(jobManagerTopic.getJMSTopic());

				// if no other jobs in 'waiting queue' (state WAITING), send to comps 
				if (jobsDb.getWaitingJobs().size() <= 1) { // that one is this job
					// forward to comp, if this fails, job is left waiting, which is ok
					compTopic.sendMessage(msg);
				}
			} 
			
			// job was not added to db, inform client
			else {
				ResultMessage resultMessage = new ResultMessage();
				resultMessage.setJobId(msg.getJobId());
				resultMessage.setState(JobState.ERROR);
				resultMessage.setErrorMessage("Could not submit job");
				endpoint.sendMessageToClientReplyChannel(clientReplyTo, resultMessage);
			}
		}

		private void handleClientCommandMessage(CommandMessage msg) throws JMSException {

			if (CommandMessage.COMMAND_GET_JOB.equals(msg.getCommand())) {
				
				Destination newClientReplyTo = msg.getReplyTo();
				String jobId = msg.getNamedParameter(ParameterMessage.PARAMETER_JOB_ID);
				
				Job job = jobsDb.updateJobReplyTo(jobId, newClientReplyTo);
				
				
				ResultMessage resultMessage = null;
				if (job != null) {
					
					// has finished with results
					if (job.getFinished() != null && job.getResults() != null) {
						resultMessage = job.getResults();
					} 
					
					// has finished, no results, for example cancelled, maybe something else too 
					else if (job.getFinished() != null && job.getResults() == null) {
						resultMessage = new ResultMessage();
						resultMessage.setJobId(jobId);
						resultMessage.setState(job.getState());
					}

					// still running
					else {
						resultMessage = new ResultMessage();
						resultMessage.setJobId(jobId);
						resultMessage.setState(JobState.RUNNING);
					}

				// job not found	
				} else {
					resultMessage = new ResultMessage();
					resultMessage.setJobId(jobId);
					resultMessage.setState(JobState.ERROR);
					resultMessage.setErrorMessage("job not found");
					logger.warn("get job for non-existent job " + jobId);
				}				
				
				endpoint.sendMessageToClientReplyChannel(newClientReplyTo, resultMessage);
				
			} else if (CommandMessage.COMMAND_CANCEL.equals(msg.getCommand())) {
				String jobId = msg.getNamedParameter(ParameterMessage.PARAMETER_JOB_ID);
				if (jobsDb.updateJobCancelled(jobsDb.getJob(jobId))) {
					compTopic.sendMessage(msg);
				}
				
			} else {
				compTopic.sendMessage(msg);
			}
		}
	}
	

	private class CompMessageListener implements MessagingListener {

		@Override
		public void onChipsterMessage(ChipsterMessage msg) {
			try {
				
				if (msg instanceof ResultMessage) {					
					handleResultMessage((ResultMessage) msg);
				} 

				else if (msg instanceof JobLogMessage) {
					handleJobLogMessage((JobLogMessage) msg);
				}

				else if (msg instanceof CommandMessage) {
					handleCompCommandMessage((CommandMessage) msg);
				}

				else {
					// TODO if has jobId, try to get client replyTo and forward there?
					logger.warn("got unexcepted message from comp: " + msg.toString());
				}
			} catch (Exception e) {
				logger.error(Exceptions.getStackTrace(e));
			}
		}

		private void handleCompCommandMessage(CommandMessage msg) throws JMSException {
			
			
			if (CommandMessage.COMMAND_OFFER.equals(msg.getCommand())) {
				String jobId = msg.getNamedParameter(ParameterMessage.PARAMETER_JOB_ID);
				String compId = msg.getNamedParameter(ParameterMessage.PARAMETER_AS_ID);
				String compHost = msg.getNamedParameter(ParameterMessage.PARAMETER_HOST);
				if (jobId == null || jobId.isEmpty() || compId == null || compId.isEmpty()) {
					logger.warn(String.format("invalid offer, jobId: %s, compId: %s", jobId, compId));
					return;
				}

				Job job = jobsDb.getJob(jobId);
				if (job == null) {
					logger.warn("offer for non-existent job " + jobId);
					return;
				}
				
				// decide whether to schedule
				boolean scheduleJob = false;

				if (job.getState() == JobState.WAITING) {
					scheduleJob = true;
				} 
				
				
//				// job has never been reported by any analysis server
//				else if (job.getSecondsSinceCreated() > JOB_DEAD_AFTER && job.getSeen() == null) { 
//					scheduleJob = true;
//				} 
//
//				// the job has not recently been reported by any analysis server
////				else if (job.getSecondsSinceLastSeen() > JOB_DEAD_AFTER) {
////					scheduleJob = true;
////				}

				// schedule
				if (scheduleJob) {
					logger.info("scheduling job " + jobId);

					// create accept message
					CommandMessage acceptMessage = new CommandMessage(CommandMessage.COMMAND_ACCEPT_OFFER);
					
					acceptMessage.addNamedParameter(ParameterMessage.PARAMETER_JOB_ID, jobId);
					acceptMessage.addNamedParameter(ParameterMessage.PARAMETER_AS_ID, compId);
					acceptMessage.setUsername(job.getJobMessage().getUsername());
					
					// TODO this needed for what?
					acceptMessage.setReplyTo(jobManagerTopic.getJMSTopic());

					// send accept
					compTopic.sendMessage(acceptMessage);

					// update job state
					jobsDb.updateJobScheduled(job, compId, compHost);
				}
				
			} else if (CommandMessage.COMMAND_COMP_AVAILABLE.equals(msg.getCommand())) {
				scheduleWaitingJobs();
			}
				
			else {
				// get job reply-to and send there
				try {
					String jobId = msg.getNamedParameter(ParameterMessage.PARAMETER_JOB_ID);
					Destination destination = jobsDb.getJob(jobId).getReplyTo();
					logger.warn("sending command message to original replyTo: " + destination);
					endpoint.sendMessageToClientReplyChannel(destination, msg);
				} catch (Exception e) {
					logger.warn("failed to forward command message to client " + msg.getCommand());
				}
			}
		}

		private void handleJobLogMessage(JobLogMessage msg) {
			jobsDb.updateJobRunning(jobsDb.getJob(msg.getJobId()));
		}

		private void handleResultMessage(ResultMessage msg) throws JMSException {
			String jobId = msg.getJobId();
			Job job = jobsDb.getJob(jobId);
			
			if (job == null) {
				logger.warn("no job found for result message with job id: " + jobId);
				return;
			}
			
			// check if already finished here in jobmanager
			if (job.getFinished() != null) {
				logger.warn(String.format("result message for already finished job %s, job state in jobmanager is %s, result message state is %s, comp id is ", jobId, job.getState(), msg.getState(), msg.getNamedParameter(ParameterMessage.PARAMETER_AS_ID)));
				return;
			}
			
			JobState jobStateFromComp = msg.getState();
			
			if (jobStateFromComp == JobState.COMPLETED ||
					jobStateFromComp == JobState.FAILED ||
					jobStateFromComp == JobState.FAILED_USER_ERROR ||
					jobStateFromComp == JobState.ERROR ||
					jobStateFromComp == JobState.TIMEOUT ||
					jobStateFromComp == JobState.CANCELLED) {
				logger.info("job " + jobId + " finished with end state: " + jobStateFromComp);
				
				// don't continue if update fails
				if (!jobsDb.updateJobFinished(job, jobStateFromComp, msg)) {
					return;
				}
				
			} else if (jobStateFromComp == JobState.RUNNING) {
				// don't continue if update fails
				if (!jobsDb.updateJobRunning(job)) {
					return;
				}
					
			} else if (jobStateFromComp == JobState.NEW) {
				// not used at the moment, no need to forward to the client
				return;
			
			} else if (jobStateFromComp == JobState.COMP_BUSY) {
				// TODO refactor away
			    // jobsDb.updateJobWaiting(jobId);
				return;
				
			} else {
				logger.warn("job " + jobId + " in state " + jobStateFromComp + ", sending result message to " + job.getReplyTo());
			}
			// if things are ok, send the result message also to the client
			endpoint.sendMessageToClientReplyChannel(job.getReplyTo(), msg);
			
		}
	}


	public JobManager(String configURL) throws Exception {
		
		// initialise dir, config and logging
		DirectoryLayout.initialiseServerLayout(
		        Arrays.asList(new String[] {"jobmanager"}), configURL);
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();

		logger = Logger.getLogger(JobManager.class);
		logger.info("starting jobmanager service...");

		jobMaxWaitTime = configuration.getInt("jobmanager", "job-max-wait-time");
		
		// initialize jobs db
		this.jobsDb = new JobManagerDB(configuration);
		
		// initialize communications
		this.endpoint = new JMSMessagingEndpoint(this);

		fromClientTopic = endpoint.createTopic(Topics.Name.AUTHORISED_REQUEST_TOPIC, AccessMode.READ);
		fromClientTopic.setListener(new ClientMessageListener());

		jobManagerAdminTopic = endpoint.createTopic(Topics.Name.JOBMANAGER_ADMIN_TOPIC, AccessMode.READ);
		jobManagerAdminTopic.setListener(this);
		
		jobManagerTopic = endpoint.createTopic(Topics.Name.JOBMANAGER_TOPIC, AccessMode.READ);
		jobManagerTopic.setListener(new CompMessageListener());
		
		compTopic = endpoint.createTopic(Topics.Name.AUTHORIZED_MANAGED_REQUEST_TOPIC, AccessMode.WRITE);

		// create keep-alive thread and register shutdown hook
		KeepAliveShutdownHandler.init(this);
		
		logger.info("jobmanager is up and running [" + ApplicationConstants.VERSION + "]");
		logger.info("[mem: " + SystemMonitorUtil.getMemInfo() + "]");
	}
	
	
	

	private void scheduleWaitingJobs() {
		List<Job> waitingJobs = jobsDb.getWaitingJobs();
		if (waitingJobs.size() > 0) {
			logger.info("rescheduling " + waitingJobs.size() + " waiting jobs");
		}
		
		List<String> jobsToBeExpired = new LinkedList<String>(); // avoid removing during iteration 

		// reschedule
		for (Job job: waitingJobs) {
			try {
				if (!rescheduleJob(job.getJobId())) {
					jobsToBeExpired.add(job.getJobId());
				};
			} catch (Exception e) {
				logger.warn("could not reschedule job " + job.getJobId(), e);
			}
		}

		// remove expired (and non-existent)
		for (String jobId: jobsToBeExpired) {
			try {
				jobsDb.updateJobMaxWaitTimeReached(jobId);

				// inform client
				Job job = jobsDb.getJob(jobId); 
				if ( job != null) {
					logger.warn("sending job wait expired for job " + jobId);
					ResultMessage msg = new ResultMessage();
					msg.setJobId(jobId);
					msg.setState(JobState.FAILED);
					msg.setErrorMessage("There was no computing server available to run this job, please try again later on");
					try {
						endpoint.sendMessageToClientReplyChannel(job.getReplyTo(), msg);
					} catch (Exception e) {
						// avoid unnecessary logging if client is not there
					}
				}
			} catch (Exception e) {
				logger.error(Exceptions.getStackTrace(e));
			}
		}
	}
	
	/**
	 * 
	 * @param jobId
	 * @return false if job should be removed (expired or non-existent)
	 */
	private boolean rescheduleJob(String jobId) {
		Job job = jobsDb.getJob(jobId);
		
		if (job == null) {
			logger.warn("trying to reschedule non existing job " + jobId);
			return false;
		}

		if (job.getSecondsSinceCreated() > jobMaxWaitTime) {
			logger.warn("max wait time reached for job " + jobId);
			return false;
		}
	
		try {
			
			JobMessage jobMessage = job.getJobMessage();
			// set replyTo to jobmanager
			jobMessage.setReplyTo(jobManagerTopic.getJMSTopic());
			
			compTopic.sendMessage(jobMessage);
		} catch (JMSException e) {
			logger.error("send message failed when reschedulig job " + jobId);
			return true; // job has not expired so not removing, try again later
		}
	
		return true;
	}

	
	
	@Override
	public void onChipsterMessage(ChipsterMessage msg) {
		try {

			if (msg instanceof ServerStatusMessage) {
				logger.warn("got ServerStatusMessage");
			} else if (msg instanceof CommandMessage) {
				CommandMessage commandMessage = (CommandMessage)msg;
				if (CommandMessage.COMMAND_LIST_RUNNING_JOBS.equals(commandMessage.getCommand())) {
					logger.info("got list running jobs");
					handleListRunningJobs(commandMessage);

				} else if (CommandMessage.COMMAND_CANCEL.equals(commandMessage.getCommand())) {
					
					String jobId = commandMessage.getNamedParameter(ParameterMessage.PARAMETER_JOB_ID);
					logger.info(String.format("cancel request from admin web for job %s", jobId));
					Job job = jobsDb.getJob(jobId);
					if (jobsDb.updateJobCancelled(job)) {
						compTopic.sendMessage(commandMessage);
						logger.info(String.format("sending cancel for job %s to comps", jobId));
						// inform client
						ResultMessage cancelMessage = new ResultMessage(job.getJobId(), JobState.CANCELLED, "cancelled by admin", null, null, null);
						endpoint.sendMessageToClientReplyChannel(job.getReplyTo(), cancelMessage);
					}

				} else if (CommandMessage.COMMAND_PURGE_OLD_JOBS.equals(commandMessage.getCommand())) {
					jobsDb.purgeOldJobs();
				
				} else if (CommandMessage.COMMAND_GET_STATUS_REPORT.equals(commandMessage.getCommand())) {
					
					CommandMessage reply = new CommandMessage();
					String sysStats = SystemMonitorUtil.getMemInfo();

					String report = "";
					report += "JOBS\n\n" +
							"waiting: " + jobsDb.getWaitingJobs().size() + "\n" +
							"running: " + jobsDb.getRunningJobs().size() + "\n" +
							"all: " + jobsDb.getJobCount() + "\n" +
							"\n";
					
					report += "MEMORY\n\n";
					report += sysStats + "\n";

					reply.addNamedParameter(ParameterMessage.PARAMETER_STATUS_REPORT, report);
					endpoint.replyToMessage(commandMessage, reply);

				} else {
					logger.warn("got unexpected command message: " + commandMessage.getCommand());
				}
			}

			else {
				logger.warn("unexpected message: " + msg.toString());
			}
		} catch (Exception e) {
			logger.error(Exceptions.getStackTrace(e));
		}
	}
	
	private void handleListRunningJobs(CommandMessage commandMessage) throws JMSException {
		ArrayList<HashMap<String, Object>> jobs = new ArrayList<>();
		List<Job> runningJobs = jobsDb.getRunningJobs();
		for (Job job : runningJobs) {
			JobMessage jobMessage = job.getJobMessage();
			String host = job.getCompHost();
			if (host == null || "".equals(host)) {
				host = job.getCompId();
			}
			JobLogMessage jobLogMessage = new JobLogMessage(jobMessage.getToolId(), job.getState(), null, job.getJobId(), job.getCreated(), job.getFinished(), null, null, jobMessage.getUsername(), host);
			jobs.add(jobLogMessage.toMap());
		}		
		String json = new Gson().toJson(jobs);
		JsonMessage reply = new JsonMessage(json);
		endpoint.replyToMessage(commandMessage, reply);
	}

	@Override
	public String getName() {
		return "jobmanager";
	}

	@Override
	public void shutdown() {
		logger.info("shutdown requested");

		// close messaging endpoint
		try {
			this.endpoint.close();
		} catch (JMSException e) {
			logger.error("closing messaging endpoint failed", e);
		}

		logger.info("shutting down");
	}

}
