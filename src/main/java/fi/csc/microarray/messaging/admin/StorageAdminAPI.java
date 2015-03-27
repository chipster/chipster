package fi.csc.microarray.messaging.admin;

import java.io.IOException;
import java.util.LinkedList;
import java.util.List;
import java.util.concurrent.CountDownLatch;

import javax.jms.JMSException;

import org.apache.log4j.Logger;
import org.joda.time.format.DateTimeFormatter;
import org.joda.time.format.ISODateTimeFormat;

import fi.csc.microarray.config.ConfigurationLoader.IllegalConfigurationException;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.SuccessMessageListener;
import fi.csc.microarray.messaging.TempTopicMessagingListenerBase;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.message.ChipsterMessage;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.ParameterMessage;
import fi.csc.microarray.messaging.message.SuccessMessage;
import fi.csc.microarray.util.Strings;

/**
 * This class uses JMS messages to send data queries and converts result messages to
 * Java objects. The methods wait for the results, turning asynchronous messages to 
 * blocking method calls.
 * 
 * @author klemela
 */
public class StorageAdminAPI extends ServerAdminAPI {
	
	@SuppressWarnings("unused")
	private static final Logger logger = Logger.getLogger(StorageAdminAPI.class);

	public interface StorageEntryListener {
		public void process(List<StorageEntry> entries);
	}
	
	public StorageAdminAPI(MessagingEndpoint endpoint) throws IOException, IllegalConfigurationException, MicroarrayException, JMSException {
		super(Topics.Name.FILEBROKER_ADMIN_TOPIC, endpoint);
	}
	
	public Long[] getStorageUsage() throws JMSException, InterruptedException {
		
		StorageTotalsMessageListener listener = new StorageTotalsMessageListener();
		return listener.query();		
	}

	public List<StorageEntry> listStorageUsageOfSessions(String username) throws JMSException, InterruptedException {
		
		StorageEntryMessageListener listener = new StorageEntryMessageListener();
		listener.query(getTopic(), username);
		return listener.getEntries();
	}

	public List<StorageAggregate> listStorageUsageOfUsers() throws JMSException, InterruptedException {
		
		StorageAggregateMessageListener listener = new StorageAggregateMessageListener();
		return listener.query();
	}		
	
	public void deleteRemoteSession(String sessionID) throws JMSException, MicroarrayException {
		SuccessMessageListener replyListener = new SuccessMessageListener();  
		
		CommandMessage removeRequestMessage = new CommandMessage(CommandMessage.COMMAND_REMOVE_SESSION);
		removeRequestMessage.addNamedParameter(ParameterMessage.PARAMETER_SESSION_UUID, sessionID); 
		getTopic().sendReplyableMessage(removeRequestMessage, replyListener);

		SuccessMessage reply = replyListener.waitForReply(TIMEOUT, TIMEOUT_UNIT);

		checkSuccessMessage(reply, "delete session");		
	}
	
	private class StorageTotalsMessageListener extends TempTopicMessagingListenerBase {

		private CountDownLatch latch;
		private Long usedSpace = null;
		private Long freeSpace = null;

		public Long[] query() throws JMSException, InterruptedException {

			latch = new CountDownLatch(1);

			try {
				CommandMessage request = new CommandMessage(CommandMessage.COMMAND_GET_STORAGE_USAGE_TOTALS);

				getTopic().sendReplyableMessage(request, this);
				latch.await(TIMEOUT, TIMEOUT_UNIT);

				if (usedSpace != null && freeSpace != null) {
					return new Long[] { usedSpace, freeSpace };
				} else {
					return null;
				}
			} finally {
				// close temp topic
				this.cleanUp();
			}
		}


		public void onChipsterMessage(ChipsterMessage msg) {
			ParameterMessage resultMessage = (ParameterMessage) msg;

			String sizesString = resultMessage.getNamedParameter(ParameterMessage.PARAMETER_SIZE_LIST);

			String[] sizes = sizesString.split("\t");

			try {
				usedSpace = Long.parseLong(sizes[0]);
				freeSpace = Long.parseLong(sizes[1]);
			} catch (Exception e) {
				usedSpace = 0L;
				freeSpace = Long.MAX_VALUE;
			}
			
			latch.countDown();
		}
	}
	
	
	public static class StorageEntryMessageListener extends TempTopicMessagingListenerBase {

		private List<StorageEntry> entries;
		private long quota;
		private long quotaWarning;
		private CountDownLatch latch;
		private long storageUsage;
		
		public void query(MessagingTopic topic, String username) throws JMSException, InterruptedException {
			
			try {
				latch = new CountDownLatch(1);

				CommandMessage request = new CommandMessage(CommandMessage.COMMAND_LIST_STORAGE_USAGE_OF_SESSIONS);
				if (username != null) {
					request.addNamedParameter("username", username);
				}

				topic.sendReplyableMessage(request, this);			
				latch.await(TIMEOUT, TIMEOUT_UNIT);
			} finally {
				cleanUp();
			}
		}
		
		public List<StorageEntry> getEntries() {	
			return entries;
		}
		
		public long getQuota() {
			return quota;
		}
		
		public long getQuotaWarning() {
			return quotaWarning;
		}
		
		public long getStorageUsage() {
			return storageUsage;
		}

		public void onChipsterMessage(ChipsterMessage msg) {
			ParameterMessage resultMessage = (ParameterMessage) msg;

			String usernamesString =  resultMessage.getNamedParameter(ParameterMessage.PARAMETER_USERNAME_LIST);
			String namesString = resultMessage.getNamedParameter(ParameterMessage.PARAMETER_SESSION_NAME_LIST);
			String sizesString = resultMessage.getNamedParameter(ParameterMessage.PARAMETER_SIZE_LIST);
			String datesString = resultMessage.getNamedParameter(ParameterMessage.PARAMETER_DATE_LIST);
			String idsString = resultMessage.getNamedParameter(ParameterMessage.PARAMETER_SESSION_UUID_LIST);
			String quotaString = resultMessage.getNamedParameter(ParameterMessage.PARAMETER_QUOTA);
			String quotaWarningString = resultMessage.getNamedParameter(ParameterMessage.PARAMETER_QUOTA_WARNING);
			String storageUsageString = resultMessage.getNamedParameter(ParameterMessage.PARAMETER_SIZE);
			
			String[] usernames = Strings.splitUnlessEmpty(usernamesString, "\t");
			String[] names = Strings.splitUnlessEmpty(namesString, "\t");
			String[] sizes = Strings.splitUnlessEmpty(sizesString, "\t");
			String[] dates = Strings.splitUnlessEmpty(datesString, "\t");
			String[] ids = Strings.splitUnlessEmpty(idsString, "\t");
			
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
			
			quota = Long.parseLong(quotaString);
			quotaWarning = Long.parseLong(quotaWarningString);
			storageUsage = Long.parseLong(storageUsageString);

			latch.countDown();
		}
	}

	private class StorageAggregateMessageListener extends TempTopicMessagingListenerBase {

		private CountDownLatch latch;
		private List<StorageAggregate> entries;

		public List<StorageAggregate> query() throws JMSException, InterruptedException {

			latch = new CountDownLatch(1);

			try {
				CommandMessage request = new CommandMessage(CommandMessage.COMMAND_LIST_STORAGE_USAGE_OF_USERS);

				getTopic().sendReplyableMessage(request, this);
				latch.await(TIMEOUT, TIMEOUT_UNIT);

				return entries;
			} finally {
				// close temp topic
				this.cleanUp();
			}
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
}
