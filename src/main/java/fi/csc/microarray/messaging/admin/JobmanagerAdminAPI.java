package fi.csc.microarray.messaging.admin;

import java.io.IOException;
import java.lang.reflect.Type;
import java.util.Collection;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.TimeUnit;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.JsonMessageListener;
import fi.csc.microarray.messaging.SuccessMessageListener;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.JobLogMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;
import fi.csc.microarray.messaging.message.SuccessMessage;

public class JobmanagerAdminAPI extends ServerAdminAPI {
		
	private static final Logger logger = Logger.getLogger(JobmanagerAdminAPI.class);
	
	public static interface JobsListener {
		public void statusUpdated(Collection<JobsEntry> collection);
	}

	public JobmanagerAdminAPI() throws IOException, IllegalConfigurationException, MicroarrayException, JMSException {
		super(Topics.Name.JOBMANAGER_ADMIN_TOPIC, "jobmanager-admin");
	}
	
	public HashMap<String, JobsEntry> queryRunningJobs() throws JMSException, InterruptedException, MicroarrayException {
		
		JsonMessageListener replyListener = new JsonMessageListener();

		CommandMessage request = new CommandMessage(CommandMessage.COMMAND_LIST_RUNNING_JOBS);								
		getTopic().sendReplyableMessage(request, replyListener);				

		String json = replyListener.waitForReply(TIMEOUT, TimeUnit.SECONDS);

		// parse json if available
		if (json != null) {
			HashMap<String, JobsEntry> jobs = new HashMap<>();

			// define parameterized type
			Type listType = new TypeToken<List<HashMap<String, String>>>() {}.getType();
			List<HashMap<String, String>> mapList = new Gson().fromJson(json, listType);										

			for (HashMap<String, String> map : mapList) {
				// originally comp sent this message as JobLogMessage
				JobLogMessage msg = new JobLogMessage();
				msg.fromMap(map);

				// convert to JobsEntry
				JobsEntry job = new JobsEntry();
				job.setJobId(msg.getJobId());
				String state = msg.getState().toString();
				if (msg.getStateDetail() != null) {
					state += " (" + msg.getStateDetail() + ")";
				}
				job.setStatus(state);
				job.setStartTime(msg.getStartTime());
				job.setCompHost(msg.getCompHost());
				job.setOperation(msg.getOperation());
				job.setUsername(msg.getUsername());				

				jobs.put(msg.getJobId(), job); // remove duplicates, because all comps respond same waiting jobs
			}					
			return jobs;
		} else {
			throw new MicroarrayException("no response from jobmanager about running jobs");
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
		}
	}
}