package fi.csc.microarray.analyser;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.net.InetAddress;
import java.net.URL;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.LinkedHashMap;
import java.util.Timer;
import java.util.TimerTask;
import java.util.UUID;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import fi.csc.microarray.analyser.AnalysisDescription.ParameterDescription;
import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.constants.ApplicationConstants;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingListener;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.MonitoredNodeBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.JobLogMessage;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.messaging.message.NamiMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;
import fi.csc.microarray.messaging.message.ResultMessage;
import fi.csc.microarray.messaging.message.JobMessage.ParameterSecurityPolicy;
import fi.csc.microarray.service.KeepAliveShutdownHandler;
import fi.csc.microarray.service.ShutdownCallback;
import fi.csc.microarray.util.MemUtil;

/**
 * Executes analysis jobs and handles input&output. Uses multithreading 
 * and thread pool.
 * 
 * @author Taavi Hupponen, Aleksi Kallio
 */
public class AnalyserServer extends MonitoredNodeBase implements MessagingListener, ResultCallback, ShutdownCallback {

	private static class InternalParameterSecurityPolicy implements ParameterSecurityPolicy {

		private static final int MAX_VALUE_LENGTH = 1000;
		
		public boolean isValueValid(String value, ParameterDescription parameterDescription) {
			
			// Check parameter size (DOS protection)
			if (value.length() > MAX_VALUE_LENGTH) {
				return false;
				
			} else {
				return true;
			}
		}
		
	}
	
	private static InternalParameterSecurityPolicy INTERNAL_PARAMETER_SECURITY_POLICY = new InternalParameterSecurityPolicy();
	private static AnalysisDescription SOURCECODE_FETCH_DESCRIPTION = new AnalysisDescription(null);
	
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
	private int receiveTimeout;
	private int scheduleTimeout;
	private int timeoutCheckInterval;
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
	private LinkedHashMap<String, AnalysisJob> receivedJobs = new LinkedHashMap<String, AnalysisJob>();
	private LinkedHashMap<String, AnalysisJob> scheduledJobs = new LinkedHashMap<String, AnalysisJob>();
	private LinkedHashMap<String, AnalysisJob> runningJobs = new LinkedHashMap<String, AnalysisJob>();
	Timer timeoutTimer;
	

	/**
	 * 
	 * @throws Exception 
	 */
	public AnalyserServer() throws Exception {
		
		// Initialise dir, config and logging
		DirectoryLayout.initialiseServerLayout(Arrays.asList(new String[] {"comp"}));
		Configuration configuration = DirectoryLayout.getInstance().getConfiguration();

		// Initialise static variables, so late because they need logging
		SOURCECODE_FETCH_DESCRIPTION.setName("Fetch sourcecode (system internal operation)");
		SOURCECODE_FETCH_DESCRIPTION.addParameter(new AnalysisDescription.ParameterDescription("tool id", "ID (technical name) of the tool", false));

		// Initialise instance variables
		this.receiveTimeout = configuration.getInt("comp", "receive-timeout");
		this.scheduleTimeout = configuration.getInt("comp", "schedule-timeout");
		this.timeoutCheckInterval = configuration.getInt("comp", "timeout-check-interval");
		this.sweepWorkDir= configuration.getBoolean("comp", "sweep-work-dir");
		this.maxJobs = configuration.getInt("comp", "max-jobs");		
		logger = Logger.getLogger(AnalyserServer.class);
		loggerJobs = Logger.getLogger("jobs");
		loggerStatus = Logger.getLogger("status");

		
		// initialize working directory
		logger.info("starting compute service...");
		this.workDir = DirectoryLayout.getInstance().getJobsDataDirBase(id);
		
		// initialise custom scripts dir
		DirectoryLayout.getInstance().getCustomScriptsDir();
		
		// initialize executor service
		this.executorService = Executors.newCachedThreadPool();

		// initialize analysis tools
		this.toolRepository = new ToolRepository(this.workDir);
		
			
		// initialize timeout checker
		timeoutTimer = new Timer(true);
		timeoutTimer.schedule(new TimeoutTimerTask(), timeoutCheckInterval, timeoutCheckInterval);
		
		
		// initialize communications
		this.endpoint = new MessagingEndpoint(this);
		
		MessagingTopic analyseTopic = endpoint.createTopic(Topics.Name.AUTHORISED_REQUEST_TOPIC, AccessMode.READ);
		analyseTopic.setListener(this);
		
		managerTopic = endpoint.createTopic(Topics.Name.MANAGER_TOPIC, AccessMode.WRITE);
		
		fileBroker = new FileBrokerClient(this.endpoint.createTopic(Topics.Name.AUTHORISED_URL_TOPIC, AccessMode.WRITE));
		
		// create keep-alive thread and register shutdown hook
		KeepAliveShutdownHandler.init(this);
		
		logger.info("analyser is up and running [" + ApplicationConstants.NAMI_VERSION + "]");
		logger.info("[mem: " + MemUtil.getMemInfo() + "]");
	}
	

	public String getName() {
		return "analyser";
	}


	/**
	 * Process incoming message. JobMessage for submitting a job, CommandMessage for canceling one.
	 * 
	 * Also operation descriptions and source codes are requested with a JobMessage. 
	 *  
	 */
	public void onNamiMessage(NamiMessage namiMessage) {
		
		// create job, request operation descriptions or source code for operation
		if (namiMessage instanceof JobMessage) {
			
			JobMessage jobMessage = (JobMessage)namiMessage;
			
			// return the operations descriptions
			if ("describe".equals(jobMessage.getAnalysisId())) {
				logger.info("sending all descriptions");
				sendReplyMessage(jobMessage, createDescriptionsMessage(jobMessage));
				return; 
			} 
			
			// return source code for an operation
			else if ("describe-operation".equals(jobMessage.getAnalysisId())) {
				sendReplyMessage(jobMessage, createSourceCodeMessage(jobMessage));
				return; 
			}
			
			//  job message
			else {
				receiveJob(jobMessage);
			}
		}  
		
		// command messages
		else if (namiMessage instanceof CommandMessage) {
			CommandMessage commandMessage = (CommandMessage)namiMessage;
			
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
							logger.info("Executing job " + job.analysis.getFullName() + ", "+ job.getId() + ", " + job.getInputMessage().getUsername()) ;
						} else {
							logger.warn("Got ACCEPT_OFFER for job which is not scheduled.");
						}
					}
				}
				
				// client chose some other as, forget this job
				else {
					logger.debug("Removing scheduled job " + jobId);
					synchronized(jobsLock) {
						AnalysisJob jobToBeForgotten = receivedJobs.remove(jobId);
						
						// job was in the receivedQueue
						if (jobToBeForgotten != null) {
							receivedJobs.remove(jobToBeForgotten);
						} 
						
						// job was scheduled
						else {
							scheduledJobs.remove(jobId);
							activeJobRemoved();
						}
					}
					
				}
			}
			
			else if (CommandMessage.COMMAND_CANCEL.equals(commandMessage.getCommand())) {
				String jobId = commandMessage.getParameters().get(0);
				
				AnalysisJob job;
				synchronized(jobsLock) {
					if (receivedJobs.containsKey(jobId)) {
						job = receivedJobs.remove(jobId);
					} else if (scheduledJobs.containsKey(jobId)) {
						job = scheduledJobs.remove(jobId);
					} else {
						job = runningJobs.get(jobId);
					}
				}
				
				if (job != null) {
					job.cancel();
				} 
			}
			updateStatus();
		}		
		
		// unknown message
		else {
			logger.error("unidentified message: " + namiMessage.getMessageID());
		}

	}


	public File getWorkDir() {
		return workDir;
	}


	public boolean shouldSweepWorkDir() {
		return sweepWorkDir;
	}


	public void removeRunningJob(AnalysisJob job) {
		String hostname = "";
		
		try {
			hostname = InetAddress.getLocalHost().getCanonicalHostName();
		} catch (UnknownHostException e1) {
			logger.warn("Could not get local hostname.");
			hostname = "";
		}
		
		char delimiter = ';';
		loggerJobs.info(job.getId() + delimiter + 
						job.getInputMessage().getAnalysisId().replaceAll("\"", "") + delimiter +
						job.getState() + delimiter + 
						job.getInputMessage().getUsername() + delimiter + 
						job.getExecutionStartTime().toString() + delimiter +
						job.getExecutionEndTime().toString() + delimiter +
						hostname);
		
		logger.debug("Analyser server removing job " + job.getId() + "(" + job.getState() + ")");
		synchronized(jobsLock) {
			this.runningJobs.remove(job.getId());
		}
		activeJobRemoved();
	
		// send message to manager
		sendJobLogMessage(job);
	}
	

	
	public void sendJobLogMessage(AnalysisJob job) {
		JobLogMessage jobLogMessage;
		String hostname = "";
		
		try {
			hostname = InetAddress.getLocalHost().getCanonicalHostName();
		} catch (UnknownHostException e1) {
			logger.warn("Could not get local hostname.");
			hostname = "";
		}
		
		jobLogMessage = new JobLogMessage(
				job.getInputMessage().getAnalysisId().replaceAll("\"", ""),
				job.getState(),
				job.getId(),
				job.getExecutionStartTime(),
				job.getExecutionEndTime(),
				job.getResultMessage().getErrorMessage(),
				job.getResultMessage().getOutputText(),
				job.getInputMessage().getUsername(),
				hostname);
		
		try {
			managerTopic.sendMessage(jobLogMessage);
		} catch (JMSException e) {
			logger.error("Could not send job log message.", e);
		}
	}
	
	
	/**
	 * This is the callback method for a job to send the result message. When a job is finished the thread
	 * running a job will clean up all the data files after calling this method. 
	 * 
	 * For this reason, all the data must be sent before this method returns.
	 * 
	 * 
	 */
	public void sendResultMessage(NamiMessage original, ResultMessage reply) {
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
	private void sendReplyMessage(final NamiMessage original, final NamiMessage reply) {
		new Thread(new Runnable() {
			public void run() {
				try {
					endpoint.replyToMessage(original, reply);
				} catch (JMSException e) {
					logger.error("Could not send message.", e);
				}
			}
		}).start();
	}


	private void activeJobRemoved() {
		synchronized (jobsLock) {
			if (!receivedJobs.isEmpty() && ((runningJobs.size() + scheduledJobs.size() < maxJobs))) {
				AnalysisJob job = receivedJobs.values().iterator().next();
				receivedJobs.remove(job.getId());
				scheduleJob(job);
			}
			this.updateStatus();
		}
	}
	
	private void receiveJob(JobMessage jobMessage) {

		
		
		// check that we can run the requested analysis
		AnalysisDescription description = null;
		try {
			description = toolRepository.getDescription(jobMessage.getAnalysisId());
		} catch (AnalysisException e) {
			logger.warn("Could not fetch description for " + jobMessage.getAnalysisId());
			ResultMessage resultMessage = new ResultMessage("", JobState.ERROR, "", "Could not load operation.", 
					"", jobMessage.getReplyTo());
			sendReplyMessage(jobMessage, resultMessage);
			return;
		}

		if (description == null) {
			logger.info("Analysis " + jobMessage.getAnalysisId() + " not found.");
			ResultMessage resultMessage = new ResultMessage("", JobState.ERROR, "", "Operation not found.", 
					"", jobMessage.getReplyTo());
			sendReplyMessage(jobMessage, resultMessage);
			return;
		}
		

		// check if requested operation is supported, if not, just ignore the request
		if (!toolRepository.supports(description.getFullName())) {
			logger.debug("Analysis " + jobMessage.getAnalysisId() + " ( " + description.getSourceResourceName() + " ) not supported, ignoring request.");
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
			
			// run later
			else {
				
				receivedJobs.put(job.getId(), job);
				// try to send the ack message
				try {
					sendAckMessage(job);
				} catch (Exception e) {
					receivedJobs.remove(job.getId());
					logger.error("Could not send ACK for job " + job.getId());
				}
			}
		}
		updateStatus();
	}

	private void scheduleJob(AnalysisJob job) {
		synchronized(jobsLock) {
			job.setScheduleTime(new Date());
			scheduledJobs.put(job.getId(), job);
		}	

		try {
			sendOfferMessage(job);
		} catch (Exception e) {
			synchronized(jobsLock) {
				scheduledJobs.remove(job.getId());
			}
			logger.error("Could not send OFFER for job " + job.getId());
		}
		updateStatus();
	}

	
	
	
	private void sendAckMessage(AnalysisJob job) throws JMSException {
		// create ack message
		CommandMessage offerMessage = new CommandMessage(CommandMessage.COMMAND_ACK);
		offerMessage.addNamedParameter(ParameterMessage.PARAMETER_AS_ID, this.id);
		offerMessage.addNamedParameter(ParameterMessage.PARAMETER_JOB_ID, job.getId());

		// try to send the message
		sendReplyMessage(job.getInputMessage(), offerMessage);

	}
	
	private void sendOfferMessage(AnalysisJob job) throws JMSException {
		// create offer message
		CommandMessage offerMessage = new CommandMessage(CommandMessage.COMMAND_OFFER);
		offerMessage.addNamedParameter(ParameterMessage.PARAMETER_AS_ID, this.id);
		offerMessage.addNamedParameter(ParameterMessage.PARAMETER_JOB_ID, job.getId());

		// try to send the message
		sendReplyMessage(job.getInputMessage(), offerMessage);
	}
	
	
	private ResultMessage createDescriptionsMessage(JobMessage requestMessage) {
		ResultMessage resultMessage = new ResultMessage("", JobState.COMPLETED, "", "", 
				"", requestMessage.getReplyTo());
		try {
			String description = toolRepository.serialiseAsStringBuffer().toString();
			URL url = fileBroker.addFile(new ByteArrayInputStream(description.getBytes()), null);
			resultMessage.addPayload(DESCRIPTION_OUTPUT_NAME, url);
		} catch (Exception e) {
			logger.error("Could not send analysis descriptions", e);
			resultMessage.setState(JobState.ERROR);
			resultMessage.setErrorMessage("Could not send analysis descriptions.");
		}
		return resultMessage;
	}


	private ResultMessage createSourceCodeMessage(JobMessage requestMessage) {
		ResultMessage resultMessage = new ResultMessage("", JobState.COMPLETED, "", "", 
				"", requestMessage.getReplyTo());
		try {
			String name = new String(requestMessage.getParameters(INTERNAL_PARAMETER_SECURITY_POLICY, SOURCECODE_FETCH_DESCRIPTION).get(0));
			logger.info("sending source code for " + name);
			String sourceCode = toolRepository.getDescription(name).getSourceCode();
			
			byte[] bytes = sourceCode.getBytes();
			if (bytes.length == 0) {
				bytes = "<empty source code>".getBytes(); // zero length bytes content would hang upload
			}
			
			URL url = fileBroker.addFile(new ByteArrayInputStream(bytes), null);
			resultMessage.addPayload(SOURCECODE_OUTPUT_NAME, url);
			
		} catch (Exception e) {
			logger.error("Could not send analysis source code", e);
			resultMessage.setState(JobState.ERROR);
			resultMessage.setErrorMessage("Could not send analysis source code.");
		}
		return resultMessage;
	}
	
	private void updateStatus() {
		synchronized(jobsLock) {
			loggerStatus.info("received jobs: " + receivedJobs.size() + 
					", scheduled jobs: " + scheduledJobs.size() + 
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

				// get old received jobs
				for (AnalysisJob job: receivedJobs.values()) {
					if ((System.currentTimeMillis() - receiveTimeout * 1000) > job.getReceiveTime().getTime()) {
						jobsToBeRemoved.add(job);
					} else {
						break;
					}
				}
				
				// remove old received jobs
				for (AnalysisJob job: jobsToBeRemoved) {
					receivedJobs.remove(job.getId());
					logger.debug("Removing old received job: " + job.getId());
					logger.debug("Jobs received: " + receivedJobs.size() + ", scheduled: " + scheduledJobs.size() + ", running: " + runningJobs.size());
				}
				
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
					logger.debug("Jobs received: " + receivedJobs.size() + ", scheduled: " + scheduledJobs.size() + ", running: " + runningJobs.size());
				}
			}
		}
	}

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
