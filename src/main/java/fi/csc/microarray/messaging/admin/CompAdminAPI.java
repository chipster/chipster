package fi.csc.microarray.messaging.admin;

import java.io.IOException;
import java.util.Collection;
import java.util.HashMap;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.SuccessMessageListener;
import fi.csc.microarray.messaging.TempTopicMessagingListenerBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.JobLogMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;
import fi.csc.microarray.messaging.message.SuccessMessage;

/**
 * This class uses JMS messages to send data queries and converts result messages to
 * Java objects. The methods wait for the results, turning asynchronous messages to 
 * blocking method calls.
 * 
 * @author klemela
 */
public class CompAdminAPI extends ServerAdminAPI {
		
	private static final Logger logger = Logger.getLogger(CompAdminAPI.class);
	
	public static interface JobsListener {
		public void statusUpdated(Collection<JobsEntry> collection);
	}

	public CompAdminAPI() throws IOException, IllegalConfigurationException, MicroarrayException, JMSException {
		super(Topics.Name.COMP_ADMIN_TOPIC, "comp-admin");
	}

	public void queryRunningJobs(JobsListener listener) throws JMSException, InterruptedException {
		new RunningJobsMessageListener(listener).query();
	}
	
	public void stopGracefullyComp(String compId) throws MicroarrayException {
		SuccessMessageListener replyListener = new SuccessMessageListener();  
				
		try {
			CommandMessage removeRequestMessage = new CommandMessage(CommandMessage.COMMAND_STOP_GRACEFULLY_COMP);
			removeRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_HOST_ID, compId); 
			getTopic().sendReplyableMessage(removeRequestMessage, replyListener);

			SuccessMessage reply = replyListener.waitForReply(TIMEOUT, TIMEOUT_UNIT);
			
			checkSuccessMessage(reply, "stop comp gracefully");
												
		} catch (JMSException e) {
			logger.error("stopping comp gracefully failed", e);
		} finally {
			replyListener.cleanUp();
		}
	}

	public void cancelJob(String jobId) throws MicroarrayException {
		SuccessMessageListener replyListener = new SuccessMessageListener();  
		
		
		try {
			CommandMessage removeRequestMessage = new CommandMessage(CommandMessage.COMMAND_CANCEL);
			removeRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_JOB_ID, jobId); 
			getTopic().sendReplyableMessage(removeRequestMessage, replyListener);

			SuccessMessage reply = replyListener.waitForReply(TIMEOUT, TIMEOUT_UNIT);
			
			checkSuccessMessage(reply, "cancel job");
		} catch (JMSException e) {
			logger.error("cancel job failed", e);
		} finally {
			replyListener.cleanUp();
		}
	}
	
	private class RunningJobsMessageListener extends TempTopicMessagingListenerBase {

		private Lock mutex = new ReentrantLock();
		private HashMap<String, JobsEntry> jobs;
		private JobsListener listener;

		public RunningJobsMessageListener(JobsListener listener) {
			this.listener = listener;
		}


		public void query() throws JMSException, InterruptedException {

			mutex.lock();
			try {				
				CommandMessage request = new CommandMessage(CommandMessage.COMMAND_LIST_RUNNING_JOBS);
				getTopic().sendReplyableMessage(request, this);
				
				jobs = new HashMap<>();
			} finally {
				mutex.unlock();
			}
		}


		public void onChipsterMessage(ChipsterMessage msg) {
			
			mutex.lock();
			try {
				
				if (msg instanceof JobLogMessage) {
					JobLogMessage jobMsg = (JobLogMessage) msg;

					JobsEntry job = new JobsEntry();
					job.setJobId(jobMsg.getJobId());
					String state = jobMsg.getState().toString();
					if (jobMsg.getStateDetail() != null) {
						state += " (" + jobMsg.getStateDetail() + ")";
					}
					job.setStatus(state);
					job.setStartTime(jobMsg.getStartTime());
					job.setCompHost(jobMsg.getCompHost());
					job.setOperation(jobMsg.getOperation());
					job.setUsername(jobMsg.getUsername());				
					
					jobs.put(jobMsg.getJobId(), job); // remove duplicates, because all comps respond same waiting jobs
					
					listener.statusUpdated(jobs.values());
				}
				
			} finally {
				mutex.unlock();
			}
		}
	}	
}