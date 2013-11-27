package fi.csc.microarray.filebroker;

import java.util.concurrent.ConcurrentHashMap;

import javax.jms.JMSException;

import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingListener;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.TempTopicMessagingListener;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.Topics.Name;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.messaging.message.ChipsterMessage;

/**
 * Endpoint for replacing JMS communication with direct method calls when the communicating programs run in the same JVM.
 * This makes it possible avoid messaging configuration and authentication and still use the same code remotely through 
 * JMS and locally with this DirectFileBrokerEndpoint.
 * 
 * At the moment it is not possible to send message to any general MessagingListener, but only to FileServer. 
 * We want to transmit the calling MessaingEndpoint, which is not supported in the current MessagingListener interface. 
 * 
 * @author klemela
 */
public class DirectFileBrokerEndpoint implements MessagingEndpoint {
	
	private static ConcurrentHashMap<String, TempTopicMessagingListener> replyListeners = new ConcurrentHashMap<>();
	private FileServer fileServer;
	private String username;

	/**
	 * 
	 * @param fileServer
	 * @param username
	 */
	public DirectFileBrokerEndpoint(FileServer fileServer, String username) {
		this.fileServer = fileServer;
		this.username = username;
	}

	@Override
	public MessagingTopic createTopic(Name topicName, AccessMode accessMode)
			throws JMSException {
		return new DirectFileBrokerTopic(fileServer, this, username);
	}

	@Override
	public void replyToMessage(ChipsterMessage original,
			ChipsterMessage reply) throws JMSException {
		replyToMessage(original, reply, null);
	}

	@Override
	public void replyToMessage(final ChipsterMessage original,
			final ChipsterMessage reply, String replyChannel) throws JMSException {
		
		Thread thread = new Thread() {
			public void run() {
				TempTopicMessagingListener listener = replyListeners.get(original.getMessageID());
				if (listener == null) {
					throw new IllegalArgumentException("Can't send reply: " + reply + ", because there isn't reply listener for original message: " + original);
				}
				listener.onChipsterMessage(reply);
				replyListeners.remove(original);
			}
		};		
		thread.start();		
	}

	@Override
	public void close() throws JMSException {
		throw new UnsupportedOperationException("not supported by DirectFileBrokerMessagingEndpoint");		}

	@Override
	public AuthenticationRequestListener getAuthenticationListener() {
		throw new UnsupportedOperationException("not supported by DirectFileBrokerMessagingEndpoint");
	}

	@Override
	public void setAuthenticationListener(
			AuthenticationRequestListener authenticationListener) {
		throw new UnsupportedOperationException("not supported by DirectFileBrokerMessagingEndpoint");		}

	@Override
	public String getSessionID() {
		throw new UnsupportedOperationException("not supported by DirectFileBrokerMessagingEndpoint");
	}

	@Override
	public void setSessionID(String sessionID) {
		throw new UnsupportedOperationException("not supported by DirectFileBrokerMessagingEndpoint");
	}

	public void putReplyListener(ChipsterMessage message, TempTopicMessagingListener replyListener) {
		replyListeners.put(message.getMessageID(), replyListener);		
	}
	
	
	public static class DirectFileBrokerTopic extends MessagingTopic {

		private FileServer fileServer;
		private String username;
		private DirectFileBrokerEndpoint directEndpoint;

		public DirectFileBrokerTopic(FileServer fileServer, DirectFileBrokerEndpoint directEndpoint, String username) throws JMSException {
			super(null, Topics.Name.FILEBROKER_TOPIC.name(), Type.MOCK, AccessMode.WRITE, directEndpoint);
			
			this.fileServer = fileServer;
			this.username = username;
			this.directEndpoint = directEndpoint;
		}
		
		public void sendMessageWithUsername(final ChipsterMessage message) {
			
			message.setUsername(username);
			
			Thread thread = new Thread() {
				public void run() {					
					fileServer.handleMessage(getEndpoint(), message);
				}
			};
			thread.start();
		}
		
		@Override
		public void sendMessage(ChipsterMessage message) throws JMSException {
			sendMessageWithUsername(message);
		}
		
		@Override
		protected void sendReplyableMessage(final ChipsterMessage message, TempTopicMessagingListener replyListener, MessagingListener authenticationListener) throws JMSException {
						
			directEndpoint.putReplyListener(message, replyListener);						
			sendMessageWithUsername(message);
		}
			
		@Override	
		public void sendReplyableMessage(fi.csc.microarray.messaging.message.ChipsterMessage message, TempTopicMessagingListener replyListener) throws JMSException {			
			sendReplyableMessage(message, replyListener, null);			
		}
		
		@Override	
		public void setListener(MessagingListener listener) throws JMSException {
			throw new UnsupportedOperationException("not supported by DirectFileBrokerTopic");
		}
				
		@Override		
		public void removeListener() throws JMSException {
			throw new UnsupportedOperationException("not supported by DirectFileBrokerTopic");
		}
	}
}