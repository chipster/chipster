package fi.csc.microarray.analyser;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.net.InetAddress;
import java.net.URL;
import java.net.UnknownHostException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Date;
import java.util.HashSet;
import java.util.LinkedHashMap;
import java.util.Timer;
import java.util.TimerTask;
import java.util.UUID;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import javax.jms.JMSException;

import org.apache.log4j.Logger;

import fi.csc.microarray.ApplicationConstants;
import fi.csc.microarray.MicroarrayException;
import fi.csc.microarray.config.Configuration;
import fi.csc.microarray.config.DirectoryLayout;
import fi.csc.microarray.config.ConfigurationLoader.OldConfigurationFormatException;
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
import fi.csc.microarray.util.Files;
import fi.csc.microarray.util.MemUtil;

/**
 * Executes analysis jobs and handles input&output. Uses multithreading 
 * and thread pool.
 * 
 * @author Taavi Hupponen, Aleksi Kallio
 */
public class AnalyserServer extends MonitoredNodeBase implements MessagingListener, ResultCallback {

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
	private String workDirBase;
	private boolean sweepWorkDir;
	private String customScriptsDirName;
	private int maxJobs;
	
	/**
	 * Id of the analyser server instance.
	 */
	private String id = UUID.randomUUID().toString();
	
	private File workDir;
	
	
	/**
	 * Available analysis'.
	 */
	private AnalysisDescriptionRepository descriptionRepository;
	private HashSet<String> supportedOperations = new HashSet<String>();
	
	
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
	
	/**
	 * Pooling of operating system processes.
	 */
	private ProcessPool processPool;

	// synchronize with this object when accessing the job maps below
	private Object jobsLock = new Object(); 
	private LinkedHashMap<String, AnalysisJob> receivedJobs = new LinkedHashMap<String, AnalysisJob>();
	private LinkedHashMap<String, AnalysisJob> scheduledJobs = new LinkedHashMap<String, AnalysisJob>();
	private LinkedHashMap<String, AnalysisJob> runningJobs = new LinkedHashMap<String, AnalysisJob>();
	Timer timeoutTimer;
	

	/**
	 * 
	 * @throws JMSException
	 * @throws IOException if creation of working directory fails.
	 * @throws MicroarrayException
	 * @throws OldConfigurationFormatException 
	 */
	public AnalyserServer() throws JMSException, IOException, MicroarrayException, OldConfigurationFormatException {
		
		// initialise dir, config and logging
		DirectoryLayout.initialiseServerLayout();
		this.descriptionRepository = new AnalysisDescriptionRepository();
		this.receiveTimeout = Integer.parseInt(Configuration.getValue("analyser", "receive_timeout"));
		this.scheduleTimeout = Integer.parseInt(Configuration.getValue("analyser", "schedule_timeout"));
		this.timeoutCheckInterval = Integer.parseInt(Configuration.getValue("analyser", "timeout_check_interval"));
		this.workDirBase = Configuration.getValue("analyser", "work_dir");
		this.sweepWorkDir= "true".equals(Configuration.getValue("analyser", "sweep_work_dir").trim());
		this.customScriptsDirName = Configuration.getValue("analyser", "customScriptsDir");
		this.maxJobs  = Integer.parseInt(Configuration.getValue("analyser", "max_jobs"));		
		logger = Logger.getLogger(AnalyserServer.class);
		loggerJobs = Logger.getLogger("jobs");
		loggerStatus = Logger.getLogger("status");

		
		// initialize working directory
		logger.info("starting compute service...");
		if (!initWorkDir()) {
			String message  = "could not initialize working directory: " + workDirBase + File.pathSeparator + this.id;
			logger.fatal(message);
			throw new IOException(message);
		}
		
		// create custom scripts dir if not exits
		File customScripts = new File(customScriptsDirName);
		if (!customScripts.exists()) {
			if (customScripts.mkdirs()) {
				logger.debug("Created custom scripts dir: " + customScripts.toString());
			} else {
				logger.warn("Could not create custom scripts dir: " + customScripts.toString());
			}
		}
		
		// initialize executor service
		executorService = Executors.newCachedThreadPool();
		
		// initialize analysis handlers
		for (String analysisHandler : Configuration.getValues("analyser", "analysis_handlers")) {
			try {
				AnalysisHandler handler = (AnalysisHandler)Class.forName(analysisHandler).newInstance();
				descriptionRepository.addAnalysisHandler(handler);
				logger.debug("initialised handler " + analysisHandler);
			} catch (Exception e) {
				// we continue, despite some stuff can not be loaded
				logger.error(e);
			}			
		}
		
		// load analysis operation scripts/classes
		
		// load descriptions of all operations, also those not supported by this instance of AS
		// this way any AS can send the descriptions when client asks for them
		ArrayList<String> allOperations = new ArrayList<String>();
		String[] configOperations = Configuration.getValues("analyser", "operations");
		allOperations.addAll(Arrays.asList(configOperations));
		
		// load additional scripts from custom-scripts
		for (File f: Files.listFilesRecursively(customScripts)) {
			allOperations.add(f.getAbsolutePath().replace(customScripts.getAbsolutePath(), ""));
		}
		
		
		if (allOperations != null) {
			for (String operation : allOperations) {
					descriptionRepository.loadOperation(operation, false);
			}
		} else {
			logger.error("No operations found on the configuration file.");
		}

		String[] allHiddenOperations = Configuration.getValues("analyser", "hidden-operations");
		if (allHiddenOperations != null) {
			for (String operation : allHiddenOperations) {
				descriptionRepository.loadOperation(operation, true);
			}
		}

		
		// read the includedOperations and excludedOperations sections to get 
		// a list of supported operations
		
		
		String[] includedOperations; 
		String[] excludedOperationsList;
		String[] hiddenOperationsList;
		HashSet<String> excludedOperations = new HashSet<String>();
		HashSet<String> hiddenOperations = new HashSet<String>();
		
		
		// get included operations
		includedOperations = Configuration.getValues("analyser", "includeOperations");
		if (includedOperations == null) {
			logger.debug("No includeOperations section, including all operations.");
			includedOperations = allOperations.toArray(new String[allOperations.size()]);
		}
	

		// get excluded operations
		excludedOperationsList = Configuration.getValues("analyser", "excludeOperations");
		if (excludedOperationsList != null) {
			for (String value : excludedOperationsList) {
				excludedOperations.add(value);
			}
			logger.debug("Excluded operations: " +  excludedOperations.toString());
		
		} else {
			logger.debug("No excludeOperations section.");
		}
		
		// get hidden operations
		hiddenOperationsList = Configuration.getValues("analyser", "hidden-operations");
		if (hiddenOperationsList != null) {
			for (String value : hiddenOperationsList) {
				hiddenOperations.add(value);
			}
		}
		
		
		// add included if not excluded
		for (String operation : includedOperations) {
			if (!excludedOperations.contains(operation)) {
				//descriptionRepository.loadOperation(operation, false);
				supportedOperations.add(operation);
			}
		}

		// add hidden operations
		if (hiddenOperationsList != null) {
			for (String operation : hiddenOperationsList) {
				if (!excludedOperations.contains(operation)) {
					//descriptionRepository.loadOperation(operation, true);
					supportedOperations.add(operation);
				}
			}
		}
		
			
		
		// initialize process pool
		processPool = new ProcessPool(this.getWorkDir());
		
		// initialize timeout checker
		timeoutTimer = new Timer(true);
		timeoutTimer.schedule(new TimeoutTimerTask(), timeoutCheckInterval, timeoutCheckInterval);
		
		
		// initialize communications
		this.endpoint = new MessagingEndpoint(this);
		
		MessagingTopic analyseTopic = endpoint.createTopic(Topics.Name.AUTHORISED_REQUEST_TOPIC, AccessMode.READ);
		analyseTopic.setListener(this);
		
		managerTopic = endpoint.createTopic(Topics.Name.MANAGER_TOPIC, AccessMode.WRITE);
		
		fileBroker = new FileBrokerClient(this.endpoint.createTopic(Topics.Name.AUTHORISED_URL_TOPIC, AccessMode.WRITE));
		
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
							logger.info("Executing job " + job.getId());
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


	public ProcessPool getProcessPool() {
		return processPool;
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
			description = descriptionRepository.getDescription(jobMessage.getAnalysisId());
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
		if (!supportedOperations.contains(description.getSourceResourceName())) {
			logger.debug("Analysis " + jobMessage.getAnalysisId() + " ( " + description.getSourceResourceName() + " ) not supported, ignoring request.");
			return;
		}

		AnalysisJob job = description.createAnalysisJob(jobMessage, this);
		
		
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
			String description = descriptionRepository.serialiseAsStringBuffer().toString();
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
			String name = new String(requestMessage.getParameters().get(0));
			logger.info("sending source code for " + name);
			String sourceCode = descriptionRepository.getDescription(name).getSourceCode();
			
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
	
	private boolean initWorkDir() {
		workDir = new File(workDirBase, id);
		if (workDir.exists()) {
			logger.error("Working directory " + workDir + " already exists, should be unique.");
			return false;
		} else {
			return workDir.mkdirs();
		}
		
		
		/*
		// if work dir exists, sweep it
		if (workDir.exists()) {
			logger.debug("working directory exists, sweeping it.");
			
			if (!Files.delTree(workDir)) {
				logger.debug("could not delete the working directory.");
				return false;
			}
			logger.debug("working directory deleted");
		}
		
		// (re)create the work dir
		if (!workDir.mkdir()) {
			logger.debug("failed to create the working directory");
			return false;
		}
		
		return true;
		*/
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
	

}
