package fi.csc.microarray.jobmanager;

import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import javax.jms.Destination;
import javax.jms.JMSException;

import org.apache.log4j.Logger;

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
import fi.csc.microarray.messaging.message.ParameterMessage;
import fi.csc.microarray.messaging.message.ResultMessage;
import fi.csc.microarray.messaging.message.ServerStatusMessage;
import fi.csc.microarray.service.KeepAliveShutdownHandler;
import fi.csc.microarray.service.ShutdownCallback;
import fi.csc.microarray.util.Exceptions;
import fi.csc.microarray.util.SystemMonitorUtil;

public class JobManager extends MonitoredNodeBase implements MessagingListener, ShutdownCallback {

	private static long JOB_DEAD_AFTER = 60; // seconds
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
					logger.info("sending " + msg.getClass() + " directly to comps");
					msg.setReplyTo(jobManagerTopic.getJMSTopic());
					compTopic.sendMessage(msg);
				}
					
			} catch (Exception e) {
				logger.error(e);
			}
		}
	
	
		private void handleClientJobMessage(JobMessage msg) throws JMSException {
			
			Destination clientReplyTo = msg.getReplyTo();
			
			try {
				// store in db, also stores the original replyTo
				jobsDb.addJob(msg);
		
				// set replyTo to jobmanager
				msg.setReplyTo(jobManagerTopic.getJMSTopic());
				
				// forward
				compTopic.sendMessage(msg);
			} catch (Exception e) {
				jobsDb.updateJobError(msg.getJobId());
				
				ResultMessage resultMessage = new ResultMessage();
				resultMessage.setJobId(msg.getJobId());
				resultMessage.setState(JobState.ERROR);
				resultMessage.setErrorMessage("Could not submit job");
				endpoint.sendMessageToClientReplyChannel(clientReplyTo, resultMessage);
				
				logger.error("could not submit job " + msg.getJobId(), e);
			}				
		}

		private void handleClientCommandMessage(CommandMessage msg) throws JMSException {

			if (CommandMessage.COMMAND_GET_JOB.equals(msg.getCommand())) {
				
				Destination newClientReplyTo = msg.getReplyTo();
				String jobId = msg.getNamedParameter(ParameterMessage.PARAMETER_JOB_ID);
				logger.info("get job " + jobId);
				
				Job job = jobsDb.updateJobReplyTo(jobId, newClientReplyTo);

				ResultMessage resultMessage = null;
				if (job != null) {
					
					// has finished with results
					if (job.getFinished() != null && job.getResults() != null) {
						resultMessage = job.getResults();
					} 
					
					// has finished, no results, supposedly cancelled
					else if (job.getFinished() != null) {
						resultMessage = new ResultMessage();
						resultMessage.setJobId(jobId);
						resultMessage.setState(JobState.CANCELLED);
					}

					// still running
					else {
						resultMessage = new ResultMessage();
						resultMessage.setJobId(jobId);
						resultMessage.setState(JobState.RUNNING);
						resultMessage.setHeartbeat(false);
					}

				// job not found	
				} else {
					resultMessage = new ResultMessage();
					resultMessage.setJobId(jobId);
					resultMessage.setState(JobState.ERROR);
					resultMessage.setErrorMessage("job not found");
					logger.info("job " + jobId + " not found");
				}				
				
				endpoint.sendMessageToClientReplyChannel(newClientReplyTo, resultMessage);
				
			} else if (CommandMessage.COMMAND_CANCEL.equals(msg.getCommand())) {
				String jobId = msg.getNamedParameter(ParameterMessage.PARAMETER_JOB_ID);
				logger.info("cancelling job " + jobId);
				if (jobsDb.updateJobCancelled(jobId)) {
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

				else if (msg instanceof ServerStatusMessage) {
					// FIXME why these come here?
					logger.info("got server status message from comp");
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
				if (jobId == null || jobId.isEmpty()) {
					logger.error("got offer with no job id from comp");
					return;
				}

				logger.info("offer message for " + jobId);
				Job job = jobsDb.getJob(jobId);
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
//				// FIXME getSecondsSinceLastSeens thros NPE
////				else if (job.getSecondsSinceLastSeen() > JOB_DEAD_AFTER) {
////					scheduleJob = true;
////				}

				// schedule
				if (scheduleJob) {
					logger.info("scheduling job " + jobId);
					String compId = msg.getNamedParameter(ParameterMessage.PARAMETER_AS_ID);

					// update job state
					try {
						jobsDb.updateJobSubmitted(jobId, compId);
					} catch (Exception e) {
						logger.warn("could not update job comp, job: " + jobId + ", comp: " + compId);
						return;
					}
					
					// create accept message
					CommandMessage acceptMessage = new CommandMessage(CommandMessage.COMMAND_ACCEPT_OFFER);
					
					acceptMessage.addNamedParameter(ParameterMessage.PARAMETER_JOB_ID, jobId);
					acceptMessage.addNamedParameter(ParameterMessage.PARAMETER_AS_ID, compId);
					acceptMessage.setUsername(job.getJobMessage().getUsername());
					
					// TODO this needed for what?
					acceptMessage.setReplyTo(jobManagerTopic.getJMSTopic());

					// send accept
					compTopic.sendMessage(acceptMessage);
				}
				
			} else if (CommandMessage.COMMAND_COMP_AVAILABLE.equals(msg.getCommand())) {
				submitWaitingJobs();
			}
				
			else {
				// get job reply-to and send there
				String jobId = msg.getNamedParameter(ParameterMessage.PARAMETER_JOB_ID);
				Destination destination = jobsDb.getJob(jobId).getReplyTo();
				logger.info("sending command message to original replyTo: " + destination);
				endpoint.sendMessageToClientReplyChannel(destination, msg);
			}
			
		}

		private void handleJobLogMessage(JobLogMessage msg) {
			jobsDb.updateJobRunning(msg.getJobId());
		}

		private void handleResultMessage(ResultMessage msg) throws JMSException {
			String jobId = msg.getJobId();
			Job job = jobsDb.getJob(jobId);
			
			// check if already cancelled here in jobmanager
			if (job.getState() == JobState.CANCELLED) {
				return;
			}
			
			JobState jobStateFromComp = msg.getState();
			
			// TODO include all finished here?
			if (jobStateFromComp == JobState.COMPLETED ||
					jobStateFromComp == JobState.FAILED ||
					jobStateFromComp == JobState.FAILED_USER_ERROR ||
					jobStateFromComp == JobState.ERROR) {

				logger.info("results ready for job " + jobId + " " + jobStateFromComp);
				// TODO try catch this?
				jobsDb.updateJobResults(jobId, jobStateFromComp, msg);
				
			} else if (jobStateFromComp == JobState.RUNNING) {
				try {
					jobsDb.updateJobRunning(jobId);
				} catch (IllegalStateException e) {
					logger.warn("could not update running job " + jobId, e);
					return;
				}
				
			} else if (jobStateFromComp == JobState.COMP_BUSY) {
				// TODO refactor away
//				try {
//					jobsDb.updateJobWaiting(jobId);
//				} catch (IllegalStateException e) {
//					logger.warn("could not update waiting job " + jobId, e);
//					return;
//				}

			} else {
				logger.info("job " + jobId + " in state " + jobStateFromComp + ", sending result message to " + job.getReplyTo());
			}
			
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
		this.jobsDb = new JobManagerDB();
		
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
	
	
	

	private void submitWaitingJobs() {
		logger.info("rescheduling " + jobsDb.getWaitingJobs().size() + " waiting jobs");
		List<String> jobsToBeRemoved = new LinkedList<String>(); // avoid removing during iteration 

		// reschedule
		for (String jobId: jobsDb.getWaitingJobs()) {
			try {
				if (!rescheduleJob(jobId)) {
					jobsToBeRemoved.add(jobId);
				};
			} catch (Exception e) {
				logger.warn("could not reschedule job " + jobId, e);
			}
		}

		// remove expired (and non-existent)
		for (String jobId: jobsToBeRemoved) {
			jobsDb.updateJobMaxWaitTimeReached(jobId);
		
			Job job = jobsDb.getJob(jobId); 
			if ( job != null) {
				logger.info("sending result for wait expired job " + jobId + " " + job.getReplyTo());
				ResultMessage msg = new ResultMessage();
				msg.setJobId(jobId);
				msg.setState(JobState.FAILED);
				msg.setErrorMessage("There was no computing server available to run this job");
				try {
					endpoint.sendMessageToClientReplyChannel(job.getReplyTo(), msg);
				} catch (JMSException e) {
					logger.info(Exceptions.getStackTrace(e));
				}
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
			logger.info("max wait time reached for job " + jobId);
			return false;
		}
		
//      if job.retries >= MAX_JOB_RETRIES:
//      job.finished = datetime.datetime.utcnow()
//      job.results = populate_job_result_body(job_id, error_msg='maximum number of job submits exceeded, available chipster-comps cannot run the job')
//      headers = populate_headers(job.reply_to, RESULT_MESSAGE)
//      self.send_to(job.reply_to, headers, job.results)
	
	
		try {
			compTopic.sendMessage(job.getJobMessage());
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
				logger.info("got ServerStatusMessage");
			} else if (msg instanceof CommandMessage) {
				CommandMessage commandMessage = (CommandMessage)msg;
				if (CommandMessage.COMMAND_LIST_RUNNING_JOBS.equals(commandMessage.getCommand())) {
					logger.info("got list running jobs");

				} else if (CommandMessage.COMMAND_CANCEL.equals(commandMessage.getCommand())) {
					logger.info("got cancel from admin web");
					String jobId = commandMessage.getNamedParameter(ParameterMessage.PARAMETER_JOB_ID);
					logger.info("cancelling job " + jobId);
					if (jobsDb.updateJobCancelled(jobId)) {
						compTopic.sendMessage(commandMessage);
						// TODO inform client?
					}

				} else if (CommandMessage.COMMAND_PURGE_OLD_JOBS.equals(commandMessage.getCommand())) {
					logger.info("got purge old");
				} else if (CommandMessage.COMMAND_GET_STATUS_REPORT.equals(commandMessage.getCommand())) {
					logger.info("got get status report");
				} else {
					logger.warn("got unexpected command message: " + commandMessage.getCommand());
				}
			}

			else {

				logger.info("unexpected message: " + msg.toString());
			}
		} catch (Exception e) {
			logger.error(Exceptions.getStackTrace(e));
		}
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
