package fi.csc.microarray.messaging.admin;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import javax.jms.JMSException;
import javax.swing.Timer;

import org.apache.log4j.Logger;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.TempTopicMessagingListenerBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.JsonMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;
import fi.csc.microarray.messaging.message.ServerStatusMessage;
import fi.csc.microarray.messaging.message.SuccessMessage;

public class ServerAdminAPI {
	
	@SuppressWarnings("unused")
	private static final Logger logger = Logger.getLogger(ServerAdminAPI.class);
	
	public static interface StatusReportListener {
		public void statusUpdated(List<ServerStatusMessage> statuses);
	}
	
	NodeBase nodeSupport;
	
	public static class AdminNodeBase extends NodeBase {
		private String name;
		public AdminNodeBase(String name) {
			this.name = name;
		}
		public String getName() {
			return name;
		}
	};
	
	public static final long TIMEOUT = 30;
	public static final TimeUnit TIMEOUT_UNIT = TimeUnit.SECONDS;

	private MessagingTopic serverAdminTopic;

	private MessagingEndpoint messagingEndpoint;

	public ServerAdminAPI(Topics.Name topic, MessagingEndpoint endpoint) throws IOException, IllegalConfigurationException, MicroarrayException, JMSException {

		this.messagingEndpoint = endpoint;
		this.serverAdminTopic = messagingEndpoint.createTopic(topic, AccessMode.WRITE);
	}
	
	/**
	 * Get first status report
	 * 
	 * Wait until first status report arrives and return it.
	 * 
	 * @return
	 * @throws JMSException
	 * @throws InterruptedException
	 */
	public String getStatusReport() throws JMSException, InterruptedException {
		BlockingStatusReportMessageListener listener = new BlockingStatusReportMessageListener();
		return listener.query();
	}
	
	/**
	 * Get arbitrary number of status reports
	 * 
	 * It's not possible to wait for all status reports, because it is not know how many servers 
	 * there are. Therefore, a listener is called always when a new status report is received.
	 * @param wait 
	 * 
	 * @param listener
	 * @throws JMSException
	 * @throws InterruptedException
	 */
	public void getStatusReports(StatusReportListener listener, int wait) throws JMSException, InterruptedException {
		new AsyncStatusReportMessageListener(listener, wait).query();
	}
	
	private class BlockingStatusReportMessageListener extends TempTopicMessagingListenerBase {

		private CountDownLatch latch;
		private String report;

		public String query() throws JMSException, InterruptedException {

			latch = new CountDownLatch(1);

			try {
				CommandMessage request = new CommandMessage(CommandMessage.COMMAND_GET_STATUS_REPORT);

				getTopic().sendReplyableMessage(request, this);
				latch.await(TIMEOUT, TIMEOUT_UNIT);

				return report;
			} finally {
				// close temp topic
				this.cleanUp();
			}
		}


		public void onChipsterMessage(ChipsterMessage msg) {
			if (msg instanceof JsonMessage) {
				JsonMessage jsonMessage = (JsonMessage)msg;
				report = jsonMessage.getJson();
			} else if (msg instanceof ParameterMessage) {
				ParameterMessage resultMessage = (ParameterMessage) msg;
				report = resultMessage.getNamedParameter(ParameterMessage.PARAMETER_STATUS_REPORT);
			}
			latch.countDown();
		}
	}
	
	/**
	 * Listener will get a list of all reports received. Reports are ordered by their time of arrival. 
	 * 
	 * @author klemela
	 */
	public class AsyncStatusReportMessageListener extends TempTopicMessagingListenerBase {

		private Lock mutex = new ReentrantLock();
		private ArrayList<ServerStatusMessage> statuses;
		private StatusReportListener listener;
		private int cleanUpDelay;

		/**
		 * @param listener
		 * @param cleanUpDelay in seconds
		 */
		public AsyncStatusReportMessageListener(StatusReportListener listener, int cleanUpDelay) {
			this.listener = listener;
			this.cleanUpDelay = cleanUpDelay;
		}


		public void query() throws JMSException, InterruptedException {

			mutex.lock();
			try {				
				statuses = new ArrayList<ServerStatusMessage>();
				
				CommandMessage request = new CommandMessage(CommandMessage.COMMAND_GET_COMP_STATUS);
				getTopic().sendReplyableMessage(request, this);
				
			} finally {
				mutex.unlock();
			}
			
			// clean up the topic simply after a fixed delay, because
			// we don't know when all the serves have replied
			Timer cleanUpTimer = new Timer(cleanUpDelay * 1000, new ActionListener() {				
				@Override
				public void actionPerformed(ActionEvent e) {
					mutex.lock();
					try {
						cleanUp();
					} finally {
						mutex.unlock();
					}		
				}
			});				
			cleanUpTimer.setRepeats(false);
			cleanUpTimer.start();
		}


		public void onChipsterMessage(ChipsterMessage msg) {
			
			mutex.lock();
			try {								
				ServerStatusMessage status = (ServerStatusMessage) msg;
				
				statuses.add(status);
				
				listener.statusUpdated(statuses);
				
			} finally {			
				mutex.unlock();				
			}
		}
	}
	
	public MessagingTopic getTopic() {
		return this.serverAdminTopic;
	}
	
	public MessagingEndpoint getEndpoint() {
		return this.messagingEndpoint;
	}
	
	public void checkSuccessMessage(SuccessMessage reply, String context) throws MicroarrayException {
		if (reply == null ) {
			throw new MicroarrayException(context + " failed, no reply before timeout");			
		} else if (!reply.success()) {
			throw new MicroarrayException(context + " failed: " + reply.getErrorMessage() + " " + reply.getDetails() + " " + reply.getExceptionString());
		}
	}
}
