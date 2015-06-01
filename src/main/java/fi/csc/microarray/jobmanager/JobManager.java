package fi.csc.microarray.jobmanager;

import java.util.Arrays;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.messaging.JMSMessagingEndpoint;
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
import fi.csc.microarray.service.KeepAliveShutdownHandler;
import fi.csc.microarray.service.ShutdownCallback;
import fi.csc.microarray.util.SystemMonitorUtil;

public class JobManager extends MonitoredNodeBase implements MessagingListener, ShutdownCallback {

	private static Logger logger;

	private MessagingEndpoint endpoint;
	
	MessagingTopic fromClientTopic;
	MessagingTopic jobManagerAdminTopic;	
	MessagingTopic jobManagerTopic;
	MessagingTopic compTopic;
	
	
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
				// TODO
				else {	
				//	                self.send_to(TOPICS['comp_topic'], frame.headers, frame.body,
				//	                             reply_to=TOPICS['jobmanager_topic'])
				}
					
			} catch (Exception e) {
				logger.error(e);
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

				// TODO what is StatusMessage
//				else if (msg instanceof StatusMessage) {
//					// do nothing?
//				}

				// TODO
				else {
//		            try:
//		                client_topic = self.resolve_reply_to(body)
//		                logger.warning("unknown msg: %s" % frame.body)
//		                self.send_to(client_topic, headers, json.dumps(frame.body),
//		                             reply_to=headers.get('reply-to'))
//		            except ReplyToResolutionException:
//		                logger.warning("unable to resolve reply_to address")
				}

			} catch (Exception e) {
				logger.error(e);
			}
		}

		private void handleCompCommandMessage(CommandMessage msg) {
			String jobID = msg.getNamedParameter(ParameterMessage.PARAMETER_JOB_ID);
			if (jobID == null || jobID.isEmpty()) {
				logger.error("jobmanger got command message with no job id from comp");
			}
			
			if (CommandMessage.COMMAND_OFFER.equals(msg.getCommand())) {
				// TODO
//	            with self.session_scope() as session:
//	                job = get_job(session, job_id)
//	            schedule_job = False
//	            if not job.submitted:  # New job, never submitted before
//	                schedule_job = True
//	            elif job.seconds_since_created() > JOB_DEAD_AFTER and not job.seen:  # Job has never reported by any analysis server
//	                schedule_job = True
//	            elif job.seconds_since_last_seen() > JOB_DEAD_AFTER:  # The job has not recently been reported by any analysis server
//	                schedule_job = True
//
//	            if schedule_job:
//	                job_id = job.job_id
//	                as_id = msg.get('as-id')
//	                with self.session_scope() as session:
//	                    try:
//	                        update_job_comp(session, job_id, as_id)
//	                    except:
//	                        logger.exception("update_job_comp failed")
//	                body = populate_msg_body('choose', as_id, job_id)
//	                headers = populate_headers(TOPICS['comp_topic'], CMD_MESSAGE,
//	                                           session_id=job.session_id,
//	                                           reply_to=TOPICS['jobmanager_topic'])
//	                self.send_to(TOPICS['comp_topic'], headers, body=body)

				
				
			} else {
				// TODO
				// get job reply-to and send there
			}
			
		}

		private void handleJobLogMessage(JobLogMessage msg) {
			// TODO Auto-generated method stub
			
		}

		private void handleResultMessage(ResultMessage msg) {
			// TODO Auto-generated method stub
			
		}
	}

	
	
	
	
	
		
	

	public JobManager(String configURL) throws Exception {
		
		// initialise dir, config and logging
		DirectoryLayout.initialiseServerLayout(
		        Arrays.asList(new String[] {"jobmanager"}), configURL);
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();

		logger = Logger.getLogger(JobManager.class);
		logger.info("starting jobmanager service...");

		// initialize communications
		this.endpoint = new JMSMessagingEndpoint(this);
		

		fromClientTopic = endpoint.createTopic(Topics.Name.AUTHORISED_REQUEST_TOPIC, AccessMode.READ);
		fromClientTopic.setListener(new ClientMessageListener());

		jobManagerAdminTopic = endpoint.createTopic(Topics.Name.JOBMANAGER_ADMIN_TOPIC, AccessMode.READ);
		jobManagerAdminTopic.setListener(new CompMessageListener());
		
		jobManagerTopic = endpoint.createTopic(Topics.Name.JOBMANAGER_TOPIC, AccessMode.READ);
		jobManagerTopic.setListener(this);
		
		compTopic = endpoint.createTopic(Topics.Name.AUTHORIZED_MANAGED_REQUEST_TOPIC, AccessMode.WRITE);

		
		//	    'comp_admin_topic': '/topic/comp-admin-topic',
		//        TOPICS['admin_topic'],

		
		

		

		// create keep-alive thread and register shutdown hook
		KeepAliveShutdownHandler.init(this);
		
		logger.info("jobmanager is up and running [" + ApplicationConstants.VERSION + "]");
		logger.info("[mem: " + SystemMonitorUtil.getMemInfo() + "]");
	}
	
	
	private void handleClientJobMessage(JobMessage msg) throws JMSException {
		// TODO add to db
		
		// forward
		msg.setReplyTo(jobManagerTopic.getJMSTopic());
		compTopic.sendMessage(msg);
	}

	private void handleClientCommandMessage(CommandMessage msg) throws JMSException {

		if (CommandMessage.COMMAND_GET_JOB.equals(msg.getCommand())) {
			// TODO
//            client_topic = headers.get('reply-to')
//            job_id = msg.get('job-id')
//            headers = populate_headers(client_topic, RESULT_MESSAGE)
//            with self.session_scope() as session:
//                try:
//                    job = update_job_reply_to(session, job_id, client_topic)
//                    if job.finished and job.results:
//                        resp_body = job.results
//                    elif job.finished:
//                        resp_body = populate_job_result_body(job_id, exit_state='CANCELLED')
//                    else:
//                        resp_body = populate_job_running_body(job.job_id)
//                except JobNotFound:
//                    logger.exception("job %s not found" % job_id)
//                    resp_body = populate_job_result_body(job_id, error_msg='job not found')
//            self.send_to(client_topic, headers, resp_body)

			
		} else if (CommandMessage.COMMAND_CANCEL.equals(msg.getCommand())) {
			// TODO
//            job_id = msg.get('parameter0')
//            session_id = headers.get('session-id')
//            self.cancel_job(job_id, session_id)
			
		} else {
			compTopic.sendMessage(msg);
		}
	}
	
	
	
	
	@Override
	public void onChipsterMessage(ChipsterMessage msg) {
		logger.info("got message: " + msg.toString());
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
