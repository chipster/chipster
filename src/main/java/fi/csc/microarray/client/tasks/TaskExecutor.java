/*
 * Created on Feb 10, 2005
 *
 */
package fi.csc.microarray.client.tasks;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

import javax.jms.Destination;
import javax.jms.JMSException;
import javax.swing.SwingUtilities;
import javax.swing.Timer;
import javax.swing.event.SwingPropertyChangeSupport;

import org.apache.log4j.Logger;

import fi.csc.chipster.tools.ngs.LocalNGSPreprocess;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.client.operation.Operation;
import fi.csc.microarray.client.tasks.Task.State;
import fi.csc.microarray.databeans.DataBean;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.filebroker.FileBrokerException;
import fi.csc.microarray.filebroker.JMSFileBrokerClient;
import fi.csc.microarray.filebroker.NotEnoughDiskSpaceException;
import fi.csc.microarray.messaging.JobState;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.TempTopicMessagingListener;
import fi.csc.microarray.messaging.TempTopicMessagingListenerBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.JobMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;
import fi.csc.microarray.messaging.message.ResultMessage;
import fi.csc.microarray.util.IOUtils.CopyProgressListener;

/**
 * Allows easy management of local or remote tasks submitted through JMS and acts as a mediator between Swing Event Dispatch Thread and
 * other threads.
 * 
 * @author Aleksi Kallio
 * 
 */
public class TaskExecutor {
	/**
	 * Logger for this class
	 */

	private static final Logger logger = Logger.getLogger(TaskExecutor.class);
	private DataManager manager;
	private FileBrokerClient fileBroker;

	private MessagingTopic requestTopic;
	private LinkedList<Task> tasks = new LinkedList<Task>();
	private LinkedList<Task> runningTasks = new LinkedList<Task>();
	private SwingPropertyChangeSupport jobExecutorStateChangeSupport;
	private boolean eventsEnabled = false;

	private class TimeoutListener implements ActionListener {
		Task taskToMonitor;

		TimeoutListener(Task taskToMonitor) {
			this.taskToMonitor = taskToMonitor;
		}

		public void actionPerformed(ActionEvent e) {
			synchronized (taskToMonitor) {
				if (!taskToMonitor.getState().isFinished()) {
					updateTaskState(taskToMonitor, State.TIMEOUT, null, -1);
				}
			}

			removeFromRunningTasks(taskToMonitor);

			// send the cancel message
			sendCancelMessage(taskToMonitor);
		}
	};

	/**
	 * Encapsulates notification to make it passable to Event Dispatch Thread.
	 */
	private class TaskExecutorChangeNotifier implements Runnable {

		private TaskExecutor parent;

		public TaskExecutorChangeNotifier(TaskExecutor parent) {
			this.parent = parent;
		}

		public void run() {
			dispatch(new PropertyChangeEvent(parent, "runningJobCount", null, getRunningTaskCount()));
		}
	}

	private enum ResultListenerState {
		WAIT_FOR_ACK, WAIT_FOR_OFFER, WAIT_FOR_STATUS, FINISHED, TIMEOUT
	};

	/**
	 * For listening to temporary result Topics.
	 *
	 */
	private class ResultMessageListener extends TempTopicMessagingListenerBase {

		Task pendingTask; 
		// should be removed
		ResultListenerState internalState;
		String asId;

		/**
		 * @param pendingTask
		 * @param taskEventListener
		 */
		public ResultMessageListener(Task pendingTask) {
			this.pendingTask = pendingTask;
			this.internalState = ResultListenerState.WAIT_FOR_ACK;
		}

		public void onChipsterMessage(ChipsterMessage msg) {
			logger.debug("Task " + pendingTask.getId() + " got message (" + msg.getMessageID() + ") of type " + msg.getClass().getName());

			// ignore everything if we (ResultListener) are already finished
			if (internalState.equals(ResultListenerState.FINISHED)) {
				return;
			}

			// also ignore everything if task is already finished
			// this happens if task is cancelled or timeouts while we are waiting
			// for messages
			if (pendingTask.getState().isFinished()) {
				logger.debug("Task " + pendingTask.getId() + " already finished, ignoring message.");
				internalState = ResultListenerState.FINISHED;
				return;
			}

			// error message can arrive at any state (real error, not failed
			// analysis)
			if (msg instanceof ResultMessage) {
				ResultMessage resultMessage = (ResultMessage) msg;
				if (JobState.ERROR.equals(resultMessage.getState())) {
					logger.debug("Task " + pendingTask.getId() + " got result message with ERROR.");
					taskFinished(State.ERROR, resultMessage.getStateDetail(), resultMessage);
					return;
				} else if (JobState.FAILED_USER_ERROR.equals(resultMessage.getState())) {
					taskFinished(State.FAILED_USER_ERROR, resultMessage.getStateDetail(), resultMessage);
					return;
				}
			}

			// ResultListener state machine
			switch (internalState) {

			case WAIT_FOR_ACK:
				if (msg instanceof CommandMessage) {
					CommandMessage commandMessage = (CommandMessage) msg;

					// got ack message
					if (CommandMessage.COMMAND_ACK.equals(commandMessage.getCommand())) {
						logger.debug("Got ACK message.");
						internalState = ResultListenerState.WAIT_FOR_OFFER;
						// TODO set timeout

					}

					// got offer message
					else if (CommandMessage.COMMAND_OFFER.equals(commandMessage.getCommand())) {

						// store as-id
						asId = commandMessage.getNamedParameter(ParameterMessage.PARAMETER_AS_ID);
						logger.debug("Got OFFER from " + asId);

						// switch state (before sending the accept, as we
						// should be ready
						// to receive status immediately after the accept
						// has been sent)
						internalState = ResultListenerState.WAIT_FOR_STATUS;

						// send accept
						CommandMessage acceptMessage = new CommandMessage(CommandMessage.COMMAND_ACCEPT_OFFER);
						acceptMessage.addNamedParameter(ParameterMessage.PARAMETER_JOB_ID, pendingTask.getId());
						acceptMessage.addNamedParameter(ParameterMessage.PARAMETER_AS_ID, asId);
						logger.debug("Sending ACCEPT_OFFER to " + asId);

						try {
							requestTopic.sendMessage(acceptMessage);
							// TODO set timeout

						} catch (JMSException e) {
							logger.error("Could not send accept message.", e);

							// usually taskFinished would pick the error message, from
							// ResultMessage, but here we use the Exception.toString()
							pendingTask.setErrorMessage(e.toString());
							taskFinished(State.ERROR, "Sending message failed", null);
						}
					}
				}
				break;

			case WAIT_FOR_OFFER:
				if (msg instanceof CommandMessage) {
					CommandMessage commandMessage = (CommandMessage) msg;
					if (CommandMessage.COMMAND_OFFER.equals(commandMessage.getCommand())) {
						// store as-id
						asId = commandMessage.getNamedParameter(ParameterMessage.PARAMETER_AS_ID);

						// switch state (before sending the accept, as we
						// should be ready
						// to receive status immediately after the accept
						// has been sent)
						internalState = ResultListenerState.WAIT_FOR_STATUS;

						// send accept
						CommandMessage acceptMessage = new CommandMessage(CommandMessage.COMMAND_ACCEPT_OFFER);
						acceptMessage.addNamedParameter(ParameterMessage.PARAMETER_JOB_ID, pendingTask.getId());
						acceptMessage.addNamedParameter(ParameterMessage.PARAMETER_AS_ID, asId);
						try {
							requestTopic.sendMessage(acceptMessage);
							// TODO set timeout
						} catch (JMSException e) {
							logger.error("Could not send accept message.", e);
							// usually taskFinished would pick the error message, from
							// ResultMessage, but here we use the Exception.toString()
							pendingTask.setErrorMessage(e.toString());
							taskFinished(State.ERROR, "Sending message failed", null);
						}
					}
				}
				break;

			case WAIT_FOR_STATUS:

				// status message
				if (msg instanceof ResultMessage) {
					ResultMessage resultMessage = (ResultMessage) msg;
					JobState jobState = resultMessage.getState();

					switch (jobState) {

					case NEW:
						// this isn't really used at the moment
						break;
					case RUNNING:
						updateTaskState(pendingTask, State.RUNNING, resultMessage.getStateDetail(), -1);
						break;
					case COMPLETED:
						updateTaskState(pendingTask, State.TRANSFERRING_OUTPUTS, null, -1);
						try {
							extractPayloads(resultMessage);
						} catch (Exception e) {
							logger.error("Getting outputs failed", e);

							// usually taskFinished would pick the error message, from
							// ResultMessage, but here we use the Exception.toString()
							pendingTask.setErrorMessage(e.toString());
							taskFinished(State.ERROR, "Transferring outputs failed", null);
							break;
						}
						taskFinished(State.COMPLETED, null, resultMessage);
						break;
					case FAILED:
						taskFinished(State.FAILED, resultMessage.getStateDetail(), resultMessage);
						break;
					case FAILED_USER_ERROR:
						taskFinished(State.FAILED_USER_ERROR, resultMessage.getStateDetail(), resultMessage);
						break;
					case TIMEOUT:
						// Task state TIMEOUT is reserved for communications timeout
						taskFinished(State.FAILED, resultMessage.getStateDetail(), resultMessage);
						break;
					case RETRY:
						if (!pendingTask.hasBeenRetried()) {
							logger.debug("Resending job " + pendingTask.getId());

							// change the id of the task to avoid analyser
							// servers getting confused
							pendingTask.changeId();

							try {
								resendJobMessage(pendingTask, resultMessage.getReplyTo());
							} catch (Exception e) {
								logger.error("Could not resend job " + pendingTask.getId(), e);
								pendingTask.setErrorMessage(e.toString());
								taskFinished(State.ERROR, "Resending job failed", null);
							}
							pendingTask.setHasBeenRetried(true);

						} else {
							logger.error("Not resending the job message for the second time " + pendingTask.getId());

							pendingTask.setErrorMessage("Resending task failed.");
							taskFinished(State.ERROR, "Retransferring data failed", null);
						}
						break;
					}
				}
				break;
			}

			// Check if the task has been cancelled or client side timeout has occured while
			// processing the message.
			//
			// Task state can be finished without ResultListenerState being finished only
			// if someone external to ResultListenerState has changed the state of the task.
			// This happens when task is cancelled or timeout occurs. In such cases, possibly
			// created databeans are removed from the DataManager
			if (pendingTask.getState().isFinished() && internalState != ResultListenerState.FINISHED) {

				// update the state of the ResultListener
				internalState = ResultListenerState.FINISHED;

				// clean up possible output databeans
				// FIXME remove output databeans from manager!
			}
		}

		private void extractPayloads(ResultMessage resultMessage) throws JMSException, MicroarrayException, IOException {
			for (String name : resultMessage.payloadNames()) {
				logger.debug("output " + name);
				URL payloadUrl = resultMessage.getPayload(name);
				InputStream payload = fileBroker.getFile(payloadUrl); 
				DataBean bean = manager.createDataBean(name, payload);
				bean.setCacheUrl(payloadUrl);
				bean.setContentChanged(false);
				pendingTask.addOutput(name, bean);
			}
		}

		/**
		 * Utility method for doing stuff that needs to be done when task finishes.
		 * 
		 * 
		 * @param state
		 * @param stateDetail
		 * @param resultMessage
		 *            may be null
		 */
		private void taskFinished(State state, String stateDetail, ResultMessage resultMessage) {

			// cleanup temp topic
			this.cleanUp();
			
			if (resultMessage != null) {
				// possible screen output
				if (resultMessage.getOutputText() != null) {
					pendingTask.setScreenOutput(resultMessage.getOutputText());
				}
				// possible error message
				if (resultMessage.getErrorMessage() != null) {
					pendingTask.setErrorMessage(resultMessage.getErrorMessage());
				}					
				
				// source code
				pendingTask.setSourceCode(resultMessage.getSourceCode());
			}

			// end time(s)
			pendingTask.setEndTime(System.currentTimeMillis());

			// update state
			updateTaskState(pendingTask, state, stateDetail, -1);

			// update internal state
			this.internalState = ResultListenerState.FINISHED;

			// remove from running
			removeFromRunningTasks(pendingTask);
		}

	}

	public TaskExecutor(MessagingEndpoint endpoint, DataManager manager) throws JMSException {
		this.manager = manager;
		this.fileBroker = new JMSFileBrokerClient(endpoint.createTopic(Topics.Name.FILEBROKER_TOPIC, AccessMode.WRITE));
		this.requestTopic = endpoint.createTopic(Topics.Name.REQUEST_TOPIC, AccessMode.WRITE);
		this.jobExecutorStateChangeSupport = new SwingPropertyChangeSupport(this);
	}

	/**
	 * For unit testing, constructs partially incomplete object.
	 */
	protected TaskExecutor(DataManager manager) throws JMSException {
		this.manager = manager;
		this.jobExecutorStateChangeSupport = new SwingPropertyChangeSupport(this);
	}

	public Task createTask(Operation operation) {
		return new Task(operation);
	}
	

	/**
	 * Non-blocking.
	 */
	public void startExecuting(Task task) throws TaskException {
		startExecuting(task, -1);
	}

	/**
	 * Starts executing task and create ResultMessageListener to receive results. ResultMessageListener will call TaskEventListener.
	 * TaskEventListener is guaranteed to be called inside Swing/AWT Event Dispatch Thread, so there can be a considerable delay between
	 * result message receiving and notification.
	 * 
	 * @param task
	 * @param taskEventListener
	 * @param timeout
	 * @throws TaskException
	 */
	public void startExecuting(final Task task, int timeout) throws TaskException {
		logger.debug("Starting task " + task.getName());

		// ugly hack for local ngs preprocess
		if (task.getOperationID().equals("LocalNGSPreprocess.java")) {
			Runnable taskRunnable = new LocalNGSPreprocess(task);
			Session.getSession().getApplication().runBlockingTask("running " + task.getNamePrettyPrinted(), taskRunnable);
			return;
		}
		
		// log parameters
		List<String> parameters;
		try {
			parameters = task.getParameters();
			logger.debug("we have " + parameters.size() + " parameters");
			for (String parameter : parameters) {
				logger.debug("parameter: " + parameter);
			}
		} catch (MicroarrayException e1) {
			logger.error("Could not log parameters.");
		}

		// set task as running (task becomes visible in the task list)
		task.setStartTime(System.currentTimeMillis());
		addToRunningTasks(task);

		// send job message (start task) in a background thread
		new Thread(new Runnable() {
			public void run() {
				try {
					JobMessage jobMessage = new JobMessage(task.getId(), task.getOperationID(), task.getParameters());

					// handle inputs
					logger.debug("adding inputs to job message");
					updateTaskState(task, State.TRANSFERRING_INPUTS, null, -1);
					int i = 0;
					for (final String name : task.getInputNames()) {

						final int fi = i;
						CopyProgressListener progressListener = new CopyProgressListener() {

							long length = task.getInput(name).getContentLength();

							public void progress(long bytes) {
								float overall = ((float)fi) / ((float)task.getInputCount());
								float current = ((float)bytes) / ((float)length);
								float total = overall + (current / ((float)task.getInputCount()));
								updateTaskState(task, State.TRANSFERRING_INPUTS, null, Math.round(total * 100f));
							}
						};
						
						
						// transfer input contents to file broker if needed
						DataBean bean = task.getInput(name);
						try {
							bean.getLock().readLock().lock();

							// bean modified, upload
							if (bean.isContentChanged()) {
								bean.setCacheUrl(fileBroker.addFile(bean.getContentByteStream(), bean.getContentLength(), progressListener)); 
								bean.setContentChanged(false);
							} 

							// bean not modified, check cache, upload if needed
							else if (bean.getCacheUrl() != null && !fileBroker.checkFile(bean.getCacheUrl(), bean.getContentLength())){
								bean.setCacheUrl(fileBroker.addFile(bean.getContentByteStream(), bean.getContentLength(), progressListener));
							}

						} finally {
							bean.getLock().readLock().unlock();
						}

						// add the possibly new url to message
						jobMessage.addPayload(name, bean.getCacheUrl());
						
						
						logger.debug("added input " + name + " to job message.");
						i++;
					}

					updateTaskState(task, State.WAITING, null, -1);
					TempTopicMessagingListener replyListener = new ResultMessageListener(task);
					logger.debug("sending job message, jobId: " + jobMessage.getJobId());

					requestTopic.sendReplyableMessage(jobMessage, replyListener);
					
				} catch (NotEnoughDiskSpaceException nedse) {
					logger.warn("received not enough disk space when uploading input", nedse);
					updateTaskState(task, State.FAILED_USER_ERROR, "Not enough disk space", -1);
					task.setErrorMessage("There was not enough disk space in Chipster server to run the task. Please try again later.");
					removeFromRunningTasks(task);
					
				} catch (Exception e) {
					// could not send job message --> task fails
					logger.error("Could not send job message.", e);
					updateTaskState(task, State.ERROR, "Sending message failed", -1);
					removeFromRunningTasks(task);
				}
			}
		}).start();
		logger.debug("task starter thread started");

		// setup timeout checker if needed
		if (timeout != -1) {
			// we'll have to timeout this task
			Timer timer = new Timer(timeout, new TimeoutListener(task));
			timer.setRepeats(false);
			timer.start();
		}
	}

	/**
	 * Blocks until result is got. Can block infinitely, if no results are sent.
	 */
	public void execute(Task task) throws TaskException {

		startExecuting(task);

		// block until it is finished
		synchronized (runningTasks) {
			while (!task.getState().isFinished()) {
				try {
					runningTasks.wait(500);
				} catch (InterruptedException e) {
				}
			}
		}
	}

	public void kill(Task task) {
		logger.debug("TaskExecutor killing task " + task.getId());

		synchronized (task) {
			// task already finished?
			if (task.getState().isFinished()) {
				logger.debug("Task already finished, no need to cancel.");
				return;
			}
			updateTaskState(task, State.CANCELLED, null, -1);
		}

		// send the cancel message
		sendCancelMessage(task);

		removeFromRunningTasks(task);
	}

	public void killAll() {
		synchronized (runningTasks) {

			// copy of runningTasks, avoid concurrent modification by kill(Task task)
			LinkedList<Task> tasksToKill = new LinkedList<Task>();
			for (Task task : runningTasks) {
				tasksToKill.add(task);
			}
			for (Task task : tasksToKill) {
				kill(task);
			}

			runningTasks.clear();
			SwingUtilities.invokeLater(new TaskExecutorChangeNotifier(this));
		}
	}

	public Collection<Task> getTasks(boolean onlyRunning, boolean showHidden) {
		synchronized (runningTasks) {
			// select if we return only running or all
			LinkedList<Task> taskList = onlyRunning ? runningTasks : tasks;

			// if we show also hidden, we can return
			if (showHidden) {
				return taskList;
			}

			// prune away hidden tasks
			LinkedList<Task> prunedTaskList = new LinkedList<Task>();
			for (Task task : taskList) {
				if (!task.isHidden()) {
					prunedTaskList.add(task);
				}
			}

			return prunedTaskList;
		}
	}

	public int getRunningTaskCount() {
		synchronized (runningTasks) {
			Collection<Task> taskList = getTasks(true, false);
			return taskList.size();
		}
	}

	/**
	 * Adds a listener for general task execution state (how many tasks are running etc).
	 */
	public void addChangeListener(PropertyChangeListener listener) {
		jobExecutorStateChangeSupport.addPropertyChangeListener(listener);
	}

	public boolean isEventsEnabled() {
		return eventsEnabled;
	}

	public void setEventsEnabled(boolean eventsEnabled) {
		this.eventsEnabled = eventsEnabled;
	}

	/**
	 * If task is already finished, state remains unmodified.
	 * 
	 * 
	 * @param task
	 * @param state
	 * @param stateDetail
	 * @param completionPercentage 
	 */
	private void updateTaskState(Task task, State state, String stateDetail, int completionPercentage) {

		// setting the state will notify Task listeners
		synchronized (task) {

			// already finished, do not change the state
			if (task.getState().isFinished()) {
				return;
			}

			// not finished yet, change state
			else {

				task.setState(state);
				if (stateDetail != null) {
					task.setStateDetail(stateDetail);
				}
				
				task.setCompletionPercentage(completionPercentage);
				
				// notify TaskExecutor listeners
				SwingUtilities.invokeLater(new TaskExecutorChangeNotifier(this));
			}
		}
	}

	private void dispatch(PropertyChangeEvent event) {
		if (eventsEnabled) {
			jobExecutorStateChangeSupport.firePropertyChange(event);
		}
	}

	protected void addToRunningTasks(Task task) {
		synchronized (runningTasks) {
			tasks.add(task);
			runningTasks.add(task);
			runningTasks.notifyAll();
		}
		SwingUtilities.invokeLater(new TaskExecutorChangeNotifier(this));
	}

	protected void removeFromRunningTasks(Task task) {
		synchronized (runningTasks) {
			runningTasks.remove(task);
			runningTasks.notifyAll();
		}
		SwingUtilities.invokeLater(new TaskExecutorChangeNotifier(this));
	}

	private void sendCancelMessage(final Task task) {
		logger.debug("Sending cancel message for " + task.getId());

		// send message in a background thread
		new Thread(new Runnable() {
			public void run() {
				try {
					// create message
					CommandMessage commandMessage = new CommandMessage(CommandMessage.COMMAND_CANCEL);
					commandMessage.addParameter(task.getId());

					// send message
					logger.debug("Sending cancel message.");
					requestTopic.sendMessage(commandMessage);

				} catch (Exception e) {
					logger.error("Could not send cancel message for " + task.getId(), e);
				}
			}
		}).start();
		logger.debug("Message cancel thread started.");
	}

	private void resendJobMessage(Task task, Destination replyTo) throws TaskException, MicroarrayException, JMSException, IOException, FileBrokerException {

		JobMessage jobMessage = new JobMessage(task.getId(), task.getOperationID(), task.getParameters());
		for (String name : task.getInputNames()) {
			DataBean bean = task.getInput(name);
			try {
				bean.getLock().readLock().lock();
				bean.setCacheUrl(fileBroker.addFile(bean.getContentByteStream(), bean.getContentLength(), null)); // no progress listening on resends 
				bean.setContentChanged(false);
			} finally {
				bean.getLock().readLock().unlock();
			}
			
			jobMessage.addPayload(name, bean.getCacheUrl()); // no progress listening on resends
		}
		jobMessage.setReplyTo(replyTo);

		logger.debug("Retry replyTo is: " + jobMessage.getReplyTo());
		requestTopic.sendMessage(jobMessage);

	}

}
