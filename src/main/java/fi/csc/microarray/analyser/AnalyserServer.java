package fi.csc.microarray.analyser;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Timer;
import java.util.TimerTask;
import java.util.UUID;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.jms.JMSException;
import javax.xml.parsers.ParserConfigurationException;

import org.apache.log4j.Logger;
import org.xml.sax.SAXException;

import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.filebroker.JMSFileBrokerClient;
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
import fi.csc.microarray.messaging.message.ModuleDescriptionMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;
import fi.csc.microarray.messaging.message.ResultMessage;
import fi.csc.microarray.messaging.message.ServerStatusMessage;
import fi.csc.microarray.messaging.message.SourceMessage;
import fi.csc.microarray.messaging.message.SuccessMessage;
import fi.csc.microarray.service.KeepAliveShutdownHandler;
import fi.csc.microarray.service.ShutdownCallback;
import fi.csc.microarray.util.SystemMonitorUtil;

/**
 * Executes analysis jobs and handles input&output. Uses multithreading 
 * and thread pool.
 * 
 * @author Taavi Hupponen, Aleksi Kallio
 */
public class AnalyserServer extends MonitoredNodeBase implements MessagingListener, ResultCallback, ShutdownCallback {

	public static final String DESCRIPTION_OUTPUT_NAME = "description";
	public static final String SOURCECODE_OUTPUT_NAME = "sourcecode";
	
	/**
	 * Loggers.
	 */
	private static Logger logger;
	private static Logger loggerJobs;
	private static Logger loggerStatus;
	
	/**
	 * Directory for storing input and output files.
	 */
	private int scheduleTimeout;
	private int offerDelay;
	private int timeoutCheckInterval;
	private int heartbeatInterval;
	private boolean sweepWorkDir;
	private int maxJobs;
	
	/**
	 * Id of the analyser server instance.
	 */
	private String id = UUID.randomUUID().toString();
	
	private File workDir;
	
	
	/**
	 * All the analysis tools.
	 */
	private ToolRepository toolRepository;
	
	
	/**
	 * Our route to messaging fabric.
	 */
	private MessagingEndpoint endpoint;
	private MessagingTopic managerTopic;
	
	private FileBrokerClient fileBroker;
	
	/**
	 * Java utility for multithreading.
	 */
	private ExecutorService executorService;
	

	// synchronize with this object when accessing the job maps below
	private Object jobsLock = new Object(); 
	private LinkedHashMap<String, AnalysisJob> scheduledJobs = new LinkedHashMap<String, AnalysisJob>();
	private LinkedHashMap<String, AnalysisJob> runningJobs = new LinkedHashMap<String, AnalysisJob>();
	private Timer timeoutTimer;
	private Timer heartbeatTimer;
	private String localFilebrokerPath;
	private String overridingFilebrokerIp;
	
	volatile private boolean stopGracefully;
	
	/**
	 * 
	 * @throws Exception 
	 */
	public AnalyserServer(String configURL) throws Exception {
		
		// initialise dir, config and logging
		DirectoryLayout.initialiseServerLayout(
		        Arrays.asList(new String[] {"comp"}), configURL);
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();

		// Initialise instance variables
		this.scheduleTimeout = configuration.getInt("comp", "schedule-timeout");
		this.offerDelay = configuration.getInt("comp", "offer-delay");
		this.timeoutCheckInterval = configuration.getInt("comp", "timeout-check-interval");
		this.heartbeatInterval = configuration.getInt("comp", "job-heartbeat-interval");
		this.sweepWorkDir= configuration.getBoolean("comp", "sweep-work-dir");
		this.maxJobs = configuration.getInt("comp", "max-jobs");
		this.localFilebrokerPath = nullIfEmpty(configuration.getString("comp", "local-filebroker-user-data-path"));
		this.overridingFilebrokerIp = nullIfEmpty(configuration.getString("comp", "overriding-filebroker-ip"));				
		
		logger = Logger.getLogger(AnalyserServer.class);
		loggerJobs = Logger.getLogger("jobs");
		loggerStatus = Logger.getLogger("status");

		
		// initialize working directory
		logger.info("starting compute service...");
		this.workDir = DirectoryLayout.getInstance().getJobsDataDirBase(id);
		
		// initialize executor service
		this.executorService = Executors.newCachedThreadPool();

		// initialize analysis tools
		this.toolRepository = new ToolRepository(this.workDir);
					
		// initialize timeout checker
		timeoutTimer = new Timer(true);
		timeoutTimer.schedule(new TimeoutTimerTask(), timeoutCheckInterval, timeoutCheckInterval);
		
		heartbeatTimer = new Timer(true);
		heartbeatTimer.schedule(new JobHeartbeatTask(), heartbeatInterval, heartbeatInterval);
		
		// initialize communications
		this.endpoint = new JMSMessagingEndpoint(this);
		
		MessagingTopic analyseTopic = endpoint.createTopic(Topics.Name.AUTHORIZED_MANAGED_REQUEST_TOPIC, AccessMode.READ);
		analyseTopic.setListener(this);
		
		managerTopic = endpoint.createTopic(Topics.Name.JOB_LOG_TOPIC, AccessMode.WRITE);
		
		MessagingTopic filebrokerAdminTopic = endpoint.createTopic(Topics.Name.COMP_ADMIN_TOPIC, AccessMode.READ);
		filebrokerAdminTopic.setListener(new CompAdminMessageListener());
		
		fileBroker = new JMSFileBrokerClient(this.endpoint.createTopic(Topics.Name.AUTHORISED_FILEBROKER_TOPIC, AccessMode.WRITE), this.localFilebrokerPath, this.overridingFilebrokerIp);
		
		
		// create keep-alive thread and register shutdown hook
		KeepAliveShutdownHandler.init(this);
		
		logger.info("analyser is up and running [" + ApplicationConstants.VERSION + "]");
		logger.info("[mem: " + SystemMonitorUtil.getMemInfo() + "]");
	}

	private String nullIfEmpty(String value) {
		if ("".equals(value.trim())) {
			return null;
		} else {
			return value;
		}
	}


	public String getName() {
		return "analyser";
	}


	/**
	 * Process incoming message. JobMessage for submitting a job, CommandMessage for canceling one.
	 * 
	 *  
	 *  
	 */
	public void onChipsterMessage(ChipsterMessage chipsterMessage) {

		// sanity check username
		if (chipsterMessage.getUsername() == null || chipsterMessage.getUsername().equals("")) {
			logger.warn("not accepting message with null or empty username");
			try {
				ResultMessage resultMessage = new ResultMessage("", JobState.ERROR, "", "Username was null or empty.", 
						"", chipsterMessage.getReplyTo());
				sendReplyMessage(chipsterMessage, resultMessage);
			} catch (Exception e) {
				logger.warn("could not send error message for null or empty username");
			}
			return;
		}
		
		// job message
		if (chipsterMessage instanceof JobMessage) {
			JobMessage jobMessage = (JobMessage)chipsterMessage;
			receiveJob(jobMessage);
		}  
		
		// command messages
		else if (chipsterMessage instanceof CommandMessage) {
			CommandMessage commandMessage = (CommandMessage)chipsterMessage;
			
			if (CommandMessage.COMMAND_ACCEPT_OFFER.equals(commandMessage.getCommand())) {
				
				// is this AS accepted?
				String acceptedId = commandMessage.getNamedParameter(ParameterMessage.PARAMETER_AS_ID);
				String jobId = commandMessage.getNamedParameter(ParameterMessage.PARAMETER_JOB_ID);
				logger.debug("ACCEPT_OFFER for analyser: " + acceptedId + " job: " + jobId);
				
				// client chose this AS to run this job
				if (this.id.equals(acceptedId)) {
					AnalysisJob job; 
					synchronized(jobsLock) {
						// check that we have the job as scheduled
						job = scheduledJobs.get(commandMessage.getNamedParameter(ParameterMessage.PARAMETER_JOB_ID));
						if (job != null) {
							scheduledJobs.remove(jobId);
							runningJobs.put(job.getId(), job);
	
							// run the job
							executorService.execute(job);
							logger.info("Executing job " + job.analysis.getDisplayName() + "(" + job.analysis.getID() + ")" + ", "+ job.getId() + ", " + job.getInputMessage().getUsername()) ;
						} else {
							logger.warn("Got ACCEPT_OFFER for job which is not scheduled.");
						}
					}
				}
				
				// client chose some other as, forget this job if we have it as scheduled
				else {
					logger.debug("Removing scheduled job " + jobId);
					synchronized(jobsLock) {

						if (scheduledJobs.containsKey(jobId)) {
							scheduledJobs.remove(jobId);
							activeJobRemoved();
						}
					}
				}
			}
			
			// Request to send descriptions
			else if (CommandMessage.COMMAND_DESCRIBE.equals(commandMessage.getCommand())) {
				if (stopGracefully) {
					return;
				}
				
				logger.info("sending all descriptions");
	            
	            // Send descriptions for all available modules
                try {
                    List<ModuleDescriptionMessage> list;
                    list = createDescriptionsMessages(commandMessage);

                    for (ModuleDescriptionMessage msg : list) {
        	            logger.info("sending descriptions for module " + msg.getModuleName());
                        sendReplyMessage(commandMessage, msg);
                    }
                } catch (Exception e) {
                    logger.error("sending descriptions message failed", e);
                }
	            return; 
			}

			// source code request
			else if (CommandMessage.COMMAND_GET_SOURCE.equals(commandMessage.getCommand())) {
				if (stopGracefully) {
					return;
				}
	            
				logger.info("sending source code");
				SourceMessage sourceMessage = createSourceCodeMessage(commandMessage);
	            if (sourceMessage != null) {
					sendReplyMessage(commandMessage, sourceMessage);
	            }
	            return;
			}			
			
			// Request to cancel a job
			else if (CommandMessage.COMMAND_CANCEL.equals(commandMessage.getCommand())) {
				String jobId = commandMessage.getParameters().get(0);
				
				cancelJob(jobId);			
			}
			updateStatus();
		}		
		
		// unknown message
		else {
			logger.error("unidentified message: " + chipsterMessage.getMessageID());
		}

	}


	private void cancelJob(String jobId) {
		AnalysisJob job;
		synchronized(jobsLock) {
			if (scheduledJobs.containsKey(jobId)) {
				job = scheduledJobs.remove(jobId);
			} else {
				job = runningJobs.remove(jobId);
			}
		}
		
		if (job != null) {
			job.cancel();
		}
	}


	public File getWorkDir() {
		return workDir;
	}


	public boolean shouldSweepWorkDir() {
		return sweepWorkDir;
	}

	public void removeRunningJob(AnalysisJob job) {
		String hostname = getHost();	
		
		char delimiter = ';';
		try {
			loggerJobs.info(job.getId() + delimiter + job.getInputMessage().getAnalysisId().replaceAll("\"", "") + delimiter
					+ job.getState() + delimiter + job.getInputMessage().getUsername() + delimiter + job.getExecutionStartTime().toString()
					+ delimiter + job.getExecutionEndTime().toString() + delimiter + hostname);
		} catch (Exception e) {
			logger.warn("got exception when logging a job to be removed", e);
		}
		logger.debug("Analyser server removing job " + job.getId() + "(" + job.getState() + ")");
		synchronized(jobsLock) {
			this.runningJobs.remove(job.getId());
		}
		activeJobRemoved();
	
		// send message to manager
		sendJobLogMessage(job);
		
		checkStopGracefully();		
	}
	

	
	private void checkStopGracefully() {
		if (stopGracefully) {
			synchronized(jobsLock) {
				if (this.scheduledJobs.size() == 0 && this.runningJobs.size() == 0) {
					shutdown();
					System.exit(0);
				}
			}
		}
	}


	public void sendJobLogMessage(AnalysisJob job) {
		JobLogMessage jobLogMessage;		
		
		jobLogMessage = jobToMessage(job);
		
		try {
			managerTopic.sendMessage(jobLogMessage);
		} catch (JMSException e) {
			logger.error("Could not send job log message.", e);
		}
	}
	
	
	private JobLogMessage jobToMessage(AnalysisJob job) {
		
		String hostname = getHost();
		
		// current jobs in admin-web may not have starTime yet
		Date startTime = job.getExecutionStartTime();
		if (startTime == null) {
			startTime = job.getScheduleTime();
		}
		if (startTime == null) {
			startTime = job.getReceiveTime();
		}
		
		JobLogMessage jobLogMessage = new JobLogMessage(
				job.getInputMessage().getAnalysisId().replaceAll("\"", ""),
				job.getState(),
				job.getStateDetail(),
				job.getId(),
				startTime,
				job.getExecutionEndTime(),
				job.getResultMessage().getErrorMessage(),
				job.getResultMessage().getOutputText(),
				job.getInputMessage().getUsername(),
				hostname);
		
		return jobLogMessage;
	}


	/**
	 * This is the callback method for a job to send the result message. When a job is finished the thread
	 * running a job will clean up all the data files after calling this method. 
	 * 
	 * For this reason, all the data must be sent before this method returns.
	 * 
	 * 
	 */
	public void sendResultMessage(ChipsterMessage original, ResultMessage reply) {
		try {
			endpoint.replyToMessage(original, reply);
		} catch (JMSException e) {
			logger.error("Could not send ResultMessage " + reply.getMessageID());
		}
		logger.info("result message sent (" + reply.getMessageID() + " " + reply.getState() + ")");
	}

	public FileBrokerClient getFileBrokerClient() {
		return this.fileBroker;
	}
	
	
	/**
	 * Sends the message in new thread.
	 * @param original
	 * @param reply
	 */
	private void sendReplyMessage(final ChipsterMessage original, final ChipsterMessage reply) {
		new Thread(new Runnable() {
			public void run() {
				try {
					endpoint.replyToMessage(original, reply);
				} catch (JMSException e) {
					// Failing is ok, if some other comp has replied quicker and
					// the TempTopic has already been deleted
					//logger.error("Could not send message.", e);
				}
			}
		}).start();
	}


	private void activeJobRemoved() {
		this.updateStatus();
	}
	
	private void receiveJob(JobMessage jobMessage) {

		logger.info("received job request from: " + jobMessage.getUsername());
		logger.info("checking if matches guest account: " + DirectoryLayout.getInstance().getConfiguration().getString("security", "guest-username"));
		
		if (jobMessage.getUsername().equals(DirectoryLayout.getInstance().getConfiguration().getString("security", "guest-username"))) {
			ResultMessage resultMessage = new ResultMessage("", JobState.FAILED_USER_ERROR, "", "Running tools is disabled for guest users.", 
					"", jobMessage.getReplyTo());
			sendReplyMessage(jobMessage, resultMessage);
			return;
		}
		
		// don't accept new jobs when shutting down
		if (stopGracefully) {
			return;
		}
		
		
		// try to find the requested operation in tool repository
		ToolDescription description = toolRepository.getDescription(jobMessage.getAnalysisId());

		// check if this instance has the requested operation
		if (description == null) {
		    return;
		}

		// check if requested operation is supported, if not, just ignore the request
		if (!toolRepository.supports(description.getID())) {
			logger.debug("analysis " + jobMessage.getAnalysisId() + " ( " + description.getToolFile() + " ) not supported, ignoring request.");
			return;
		}

		AnalysisJob job;
		try {
			job = description.createAnalysisJob(jobMessage, this);
			
		} catch (AnalysisException e) {
			logger.warn("could not create analysis job for " + jobMessage.getAnalysisId());
			ResultMessage resultMessage = new ResultMessage("", JobState.ERROR, "", "Could not initialise operation.", 
					"", jobMessage.getReplyTo());
			sendReplyMessage(jobMessage, resultMessage);
			return;
		}
		
		// now we know that we can run this job
		// check if we could run it now or later
		synchronized(jobsLock) {
			job.setReceiveTime(new Date());
			
			
			// could run it now
			if (runningJobs.size() + scheduledJobs.size() < maxJobs) {
				scheduleJob(job);
			}
			
			// no slot to run it now, ignore it
			else {
				ResultMessage resultMessage = new ResultMessage(jobMessage.getJobId(), JobState.COMP_BUSY, "", "", "", jobMessage.getReplyTo());
				sendReplyMessage(jobMessage, resultMessage);
				return;
			}
		}
		updateStatus();
	}

	private void scheduleJob(final AnalysisJob job) {
		synchronized(jobsLock) {
			job.setScheduleTime(new Date());
			scheduledJobs.put(job.getId(), job);
		}	

		try {
			// delaying sending of the offer message can be used for
			// prioritising comp instances 
			int delay = offerDelay * (runningJobs.size() + scheduledJobs.size()-1);
			if (delay > 0 ) {
				Timer timer = new Timer("offer-delay-timer", true);
				timer.schedule(new TimerTask() {

					@Override
					public void run() {
						try {
							sendOfferMessage(job);
						} catch (JMSException e) {
							synchronized(jobsLock) {
								scheduledJobs.remove(job.getId());
							}
							logger.error("Could not send OFFER for job " + job.getId());
						}
						updateStatus();
					}

				}, delay);
			} else {
				sendOfferMessage(job);
			}
		} catch (Exception e) {
			synchronized(jobsLock) {
				scheduledJobs.remove(job.getId());
			}
			logger.error("Could not send OFFER for job " + job.getId());
		}
		updateStatus();
	}

	private void sendOfferMessage(AnalysisJob job) throws JMSException {
		// create offer message
		CommandMessage offerMessage = new CommandMessage(CommandMessage.COMMAND_OFFER);
		offerMessage.addNamedParameter(ParameterMessage.PARAMETER_AS_ID, this.id);
		offerMessage.addNamedParameter(ParameterMessage.PARAMETER_JOB_ID, job.getId());

		// try to send the message
		sendReplyMessage(job.getInputMessage(), offerMessage);
	}
	
	private List<ModuleDescriptionMessage>
	        createDescriptionsMessages(CommandMessage requestMessage)
	        throws IOException, SAXException, ParserConfigurationException {
	    List<ModuleDescriptionMessage> list = toolRepository.getModuleDescriptions();
	    for (ModuleDescriptionMessage descriptionMsg : list) {
	        descriptionMsg.setReplyTo(requestMessage.getReplyTo());
	    }
	    return list;
	}

	private SourceMessage createSourceCodeMessage(CommandMessage requestMessage) {
			String toolID = new String(requestMessage.getParameters().get(0));
			
			logger.info("sending source code for " + toolID);
			String sourceCode = toolRepository.getDescription(toolID).getSourceCode();

			if (sourceCode != null) {
				return new SourceMessage(sourceCode);
			} else {
				return null;
			}
	}
	
	private void updateStatus() {
		synchronized(jobsLock) {
			loggerStatus.info("scheduled jobs: " + scheduledJobs.size() + 
					", running jobs: " + runningJobs.size());
		}
	}
	
	/**
	 * The order of the jobs in the receivedJobs and scheduledJobs is FIFO. Because of synchronizations 
	 * this does not necessarily strictly correspond to the receiveTime and scheduleTime fields of the
	 * jobs, but is close enough.
	 * 
	 * As the jobs are ordered, it is enough to check the jobs until the first new enough job is found
	 * as the following jobs are newer (almost always, see above).
	 * 
	 * TODO send BUSY if timeout?
	 * 
	 */
	private class TimeoutTimerTask extends TimerTask {
		
		@Override
		public void run() {
			synchronized(jobsLock) {
				
				ArrayList<AnalysisJob> jobsToBeRemoved = new ArrayList<AnalysisJob>();

				// get old scheduled jobs	
				jobsToBeRemoved.clear();
				for (AnalysisJob job: scheduledJobs.values()) {
					if ((System.currentTimeMillis() - scheduleTimeout * 1000) > job.getScheduleTime().getTime()) {
						jobsToBeRemoved.add(job);
					} else {
						break;
					}
				}

				// remove old scheduled jobs
				for (AnalysisJob job: jobsToBeRemoved) {
					scheduledJobs.remove(job.getId());
					logger.debug("Removing old scheduled job: " + job.getId());
					activeJobRemoved();
				}
			}
		}
	}
	
	public class JobHeartbeatTask extends TimerTask {

		@Override
		public void run() {
			synchronized (jobsLock) {
				for (AnalysisJob job : getAllJobs()) {
					job.updateStateToClient();
				}
			}
		}	
	}
	

	/* 
	 * Service wrapper's "stop" command
	 */
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
	
	private synchronized ArrayList<AnalysisJob> getAllJobs() {
		ArrayList<AnalysisJob> allJobs = new ArrayList<AnalysisJob>();
		allJobs.addAll(scheduledJobs.values());
		allJobs.addAll(runningJobs.values());	

		return allJobs;
	}
	
	private class CompAdminMessageListener implements MessagingListener {		

		@Override
		public void onChipsterMessage(ChipsterMessage msg) {						

			try {

				if (msg instanceof CommandMessage && CommandMessage.COMMAND_GET_COMP_STATUS.equals(((CommandMessage)msg).getCommand())) {

					CommandMessage requestMessage = (CommandMessage) msg;								
																							
					ServerStatusMessage reply = SystemMonitorUtil.getSystemStats(AnalyserServer.this.workDir);

					if (stopGracefully) {
						reply.setStatus("Stopping gracefully...");
					}
					reply.setScheduledJobs(scheduledJobs.size());
					reply.setRunningJobs(runningJobs.size());
					reply.setHost(getHost());
					reply.setHostId(id);									

					endpoint.replyToMessage(requestMessage, reply);
				}

				else if (msg instanceof CommandMessage && CommandMessage.COMMAND_LIST_RUNNING_JOBS.equals(((CommandMessage)msg).getCommand())) {
					
					
					CommandMessage requestMessage = (CommandMessage) msg;
					
					synchronized (jobsLock) {
									
						for (AnalysisJob job : getAllJobs()) {
							JobLogMessage reply = jobToMessage(job);
							endpoint.replyToMessage(requestMessage, reply);
						}						
					}
				}


				else if (msg instanceof CommandMessage && CommandMessage.COMMAND_STOP_GRACEFULLY_COMP.equals(((CommandMessage)msg).getCommand())) {
					CommandMessage requestMessage = (CommandMessage) msg;
					
					String compId = ((ParameterMessage)msg).getNamedParameter(ParameterMessage.PARAMETER_HOST_ID);
					
					if (id.equals(compId)) {

						stopGracefully = true;						

						logger.info("Server received a shutdown request. "
								+ "It will shutdown after all scheduled and running jobs are completed. "
								+ "Scheduling of new jobs is disabled.");

						SuccessMessage reply = new SuccessMessage(true);									

						endpoint.replyToMessage(requestMessage, reply);
						
						checkStopGracefully();
					}
				}
				
				else if (msg instanceof CommandMessage && CommandMessage.COMMAND_CANCEL.equals(((CommandMessage)msg).getCommand())) {
					CommandMessage requestMessage = (CommandMessage) msg;
					
					String jobId = ((ParameterMessage)msg).getNamedParameter(ParameterMessage.PARAMETER_JOB_ID);
					
					logger.info("Request from CompAdminAPI to cancel a job " + jobId);
					
					cancelJob(jobId);

					SuccessMessage reply = new SuccessMessage(true);									

					endpoint.replyToMessage(requestMessage, reply);
					
					updateStatus();
				}

			} catch (Exception e) {
				logger.error(e, e);
			}
		}	
	}
}
