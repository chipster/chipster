package fi.csc.chipster.web.adminweb.data;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;

import javax.jms.JMSException;

import org.apache.log4j.Logger;
import org.joda.time.format.DateTimeFormatter;
import org.joda.time.format.ISODateTimeFormat;

import fi.csc.chipster.web.adminweb.ChipsterConfiguration;
import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.JMSMessagingEndpoint;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.SuccessMessageListener;
import fi.csc.microarray.messaging.TempTopicMessagingListenerBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;
import fi.csc.microarray.messaging.message.SuccessMessage;

/**
 * This class uses JMS messages to send data queries and converts result messages to
 * Java objects. The methods wait for the results, turning asynchronous messages to 
 * blocking method calls.
 * 
 * @author klemela
 */
public class StorageAdminAPI {
	
	private static final Logger logger = Logger.getLogger(StorageAdminAPI.class);

	public interface StorageEntryListener {
		public void process(List<StorageEntry> entries);
	}

	NodeBase nodeSupport = new NodeBase() {
		public String getName() {
			return "admin";
		}
	};
	
	private static final long TIMEOUT = 30;
	private final TimeUnit TIMEOUT_UNIT = TimeUnit.SECONDS;

	private MessagingTopic filebrokerAdminTopic;

	private MessagingEndpoint messagingEndpoint;

	public StorageAdminAPI() throws IOException, IllegalConfigurationException, MicroarrayException, JMSException {

		ChipsterConfiguration.init();
		messagingEndpoint = new JMSMessagingEndpoint(nodeSupport);
		filebrokerAdminTopic = messagingEndpoint.createTopic(Topics.Name.FILEBROKER_ADMIN_TOPIC, AccessMode.WRITE);
	}

	public List<StorageEntry> listStorageUsageOfSessions(String username) throws JMSException, InterruptedException {
		
		StorageEntryMessageListener listener = new StorageEntryMessageListener();
		return listener.query(username);
	}

	public List<StorageAggregate> listStorageUsageOfUsers() throws JMSException, InterruptedException {
		
		StorageAggregateMessageListener listener = new StorageAggregateMessageListener();
		return listener.query();
	}
	
	public void deleteRemoteSession(String sessionID) throws JMSException {
		SuccessMessageListener replyListener = new SuccessMessageListener();  
		
		
		try {
			CommandMessage removeRequestMessage = new CommandMessage(CommandMessage.COMMAND_REMOVE_SESSION);
			removeRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_SESSION_UUID, sessionID); 
			filebrokerAdminTopic.sendReplyableMessage(removeRequestMessage, replyListener);

			SuccessMessage reply = replyListener.waitForReply(TIMEOUT, TIMEOUT_UNIT);
			
			if (reply == null || !reply.success()) {
				// FIX ME, communicate properly
				throw new JMSException("failed to remove session");
			}
			
		} finally {
			replyListener.cleanUp();
		}
	}
	
	
	private class StorageEntryMessageListener extends TempTopicMessagingListenerBase {

		private List<StorageEntry> entries;
		private CountDownLatch latch;
		
		public List<StorageEntry> query(String username) throws JMSException, InterruptedException {
			
			latch = new CountDownLatch(1);
			
			CommandMessage request = new CommandMessage(CommandMessage.COMMAND_LIST_STORAGE_USAGE_OF_SESSIONS);
			request.addNamedParameter("username", username);

			filebrokerAdminTopic.sendReplyableMessage(request, this);			
			latch.await(TIMEOUT, TIMEOUT_UNIT);
			
			return entries;
		}

		public void onChipsterMessage(ChipsterMessage msg) {
			ParameterMessage resultMessage = (ParameterMessage) msg;

			String usernamesString =  resultMessage.getNamedParameter(ParameterMessage.PARAMETER_USERNAME_LIST);
			String namesString = resultMessage.getNamedParameter(ParameterMessage.PARAMETER_SESSION_NAME_LIST);
			String sizesString = resultMessage.getNamedParameter(ParameterMessage.PARAMETER_SIZE_LIST);
			String datesString = resultMessage.getNamedParameter(ParameterMessage.PARAMETER_DATE_LIST);
			String idsString = resultMessage.getNamedParameter(ParameterMessage.PARAMETER_SESSION_UUID_LIST);

			
			String[] usernames = usernamesString.split("\t");
			String[] names = namesString.split("\t");
			String[] sizes = sizesString.split("\t");
			String[] dates = datesString.split("\t");
			String[] ids = idsString.split("\t");
			
			DateTimeFormatter dateTimeFormatter = ISODateTimeFormat.dateTime();
			entries = new LinkedList<StorageEntry>();
			for (int i = 0; i < names.length; i++) {

				StorageEntry entry = new StorageEntry();
				entry.setDate(dateTimeFormatter.parseDateTime(dates[i]).toDate());
				entry.setUsername(usernames[i]);
				entry.setSize(Long.parseLong(sizes[i]));
				entry.setName(names[i]);
				entry.setID(ids[i]);
				entries.add(entry);
			}

			latch.countDown();
		}
	}	

	private class StorageAggregateMessageListener extends TempTopicMessagingListenerBase {

		private CountDownLatch latch;
		private List<StorageAggregate> entries;

		public List<StorageAggregate> query() throws JMSException, InterruptedException {

			latch = new CountDownLatch(1);

			CommandMessage request = new CommandMessage(CommandMessage.COMMAND_LIST_STORAGE_USAGE_OF_USERS);

			filebrokerAdminTopic.sendReplyableMessage(request, this);
			latch.await(TIMEOUT, TIMEOUT_UNIT);

			return entries;
		}


		public void onChipsterMessage(ChipsterMessage msg) {
			ParameterMessage resultMessage = (ParameterMessage) msg;

			String namesString = resultMessage.getNamedParameter(ParameterMessage.PARAMETER_USERNAME_LIST);
			String sizesString = resultMessage.getNamedParameter(ParameterMessage.PARAMETER_SIZE_LIST);

			String[] names = namesString.split("\t");
			String[] sizes = sizesString.split("\t");

			entries = new LinkedList<StorageAggregate>();
			for (int i = 0; i < names.length && i < sizes.length; i++) {

				StorageAggregate entry = new StorageAggregate();
				entry.setUsername(names[i]);
				entry.setSize(Long.parseLong(sizes[i]));				
				entries.add(entry);
			}

			latch.countDown();
		}
	}

	public void clean() {
		if (filebrokerAdminTopic != null) {
			try {
				filebrokerAdminTopic.delete();
			} catch (JMSException e) {
				logger.error(e);
			}
		}
		if (messagingEndpoint != null) {
			try {
				messagingEndpoint.close();
			} catch (JMSException e) {
				logger.error(e);
			}
		}
	}
}
