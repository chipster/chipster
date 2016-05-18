package fi.csc.microarray.messaging.admin;

import java.io.IOException;
import java.lang.reflect.Type;
import java.util.Collection;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.concurrent.TimeUnit;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import com.google.gson.Gson;
import com.google.gson.reflect.TypeToken;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.AuthCancelledException;
import fi.csc.microarray.messaging.JsonMessageListener;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.JobLogMessage;
import fi.csc.microarray.messaging.message.JsonMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;

public class JobmanagerAdminAPI extends ServerAdminAPI {
		
	private static final Logger logger = Logger.getLogger(JobmanagerAdminAPI.class);
	
	public static interface JobsListener {
		public void statusUpdated(Collection<JobsEntry> collection);
	}

	public JobmanagerAdminAPI(MessagingEndpoint endpoint) throws IOException, IllegalConfigurationException, MicroarrayException, JMSException {
		super(Topics.Name.JOBMANAGER_ADMIN_TOPIC, endpoint);
	}
	
	public HashMap<String, JobsEntry> queryRunningJobs() throws JMSException, InterruptedException, MicroarrayException, AuthCancelledException {
		
		JsonMessageListener replyListener = new JsonMessageListener();

		CommandMessage request = new CommandMessage(CommandMessage.COMMAND_LIST_RUNNING_JOBS);								
		getTopic().sendReplyableMessage(request, replyListener);				

		JsonMessage jsonMessage = replyListener.waitForReply(TIMEOUT, TimeUnit.SECONDS);

		// parse json if available
		if (jsonMessage != null) {
			String json = jsonMessage.getJson();
			HashMap<String, JobsEntry> jobs = new HashMap<>();

			// define parameterized type
			Type listType = new TypeToken<List<HashMap<String, String>>>() {}.getType();
			List<HashMap<String, String>> mapList = new Gson().fromJson(json, listType);										

			for (HashMap<String, String> map : mapList) {
				// originally comp sent this message as JobLogMessage
				JobLogMessage msg = new JobLogMessage();
				
				// remove keys with empty values (empty string can't be parsed to date) 
				Iterator<String> iter = map.keySet().iterator();
				while (iter.hasNext()) {
					String key = iter.next();
					if (map.get(key) == null || map.get(key).isEmpty()) {
						iter.remove();
					}
				}
				
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
				
		try {
			CommandMessage removeRequestMessage = new CommandMessage(CommandMessage.COMMAND_CANCEL);
			removeRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_JOB_ID, jobId); 
			getTopic().sendMessage(removeRequestMessage);
		} catch (JMSException e) {
			logger.error("cancel job failed", e);
		}
	}

	public void purge() {
		try {
			CommandMessage request = new CommandMessage(CommandMessage.COMMAND_PURGE_OLD_JOBS); 
			getTopic().sendMessage(request);
		} catch (JMSException e) {
			logger.error("purge old jobs failed", e);
		}
	}
}