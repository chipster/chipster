package fi.csc.microarray.messaging;

import java.util.HashMap;

import javax.jms.Destination;
import javax.jms.JMSException;

import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.Topics.Name;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.messaging.message.ChipsterMessage;

/**
 * Endpoint for replacing JMS communication with direct method calls when the communicating programs run in the same JVM.
 * This makes it possible avoid messaging configuration and authentication and still use the same code remotely through 
 * JMS and locally with this DirectMessagingEndpoint.
 * 
 * If the receiver needs to communicate simultaneously through multiple endpoints, it must implement
 * DirectMessagingListener instead of MessagingListener to know where to reply.
 * 
 * @author klemela
 */

public class DirectMessagingEndpoint implements MessagingEndpoint {
	
	private String username;

	public DirectMessagingEndpoint() {		
	}
	
	/**
	 * @param username sender of the message
	 */
	public DirectMessagingEndpoint(String username) {
		this.username = username;
	}

	private HashMap<ChipsterMessage, DirectMessagingTopic> replyTopicMap = new HashMap<ChipsterMessage, DirectMessagingTopic>();
	private HashMap<Name, DirectMessagingTopic> topicMap = new HashMap<Name, DirectMessagingTopic>();
	
	public void mapReplyTopic(ChipsterMessage original, DirectMessagingTopic replyTopic) {
		replyTopicMap.put(original, replyTopic);
	}
	
	@Override
	public MessagingTopic createTopic(Name topicName, AccessMode accessMode) throws JMSException {

		// replace unauthenticated with authencated, effectively bypassing authentication
		if (topicName == Name.FILEBROKER_TOPIC) {
			topicName = Name.AUTHORISED_FILEBROKER_TOPIC;
		} else if (topicName == Name.FEEDBACK_TOPIC) {
			topicName = Name.AUTHORISED_FEEDBACK_TOPIC;
		} else if (topicName == Name.REQUEST_TOPIC) {
			topicName = Name.AUTHORISED_REQUEST_TOPIC;
		}
		
		// add to map, if does not exist
		if (!topicMap.containsKey(topicName)) {
			topicMap.put(topicName, new DirectMessagingTopic(topicName.name(), accessMode, this));
		}
		
		return topicMap.get(topicName);
	}

	@Override
	public void replyToMessage(ChipsterMessage original, ChipsterMessage reply)	throws JMSException {
		DirectMessagingTopic replyTopic = replyTopicMap.get(original);
		if (replyTopic == null) {
			throw new IllegalArgumentException("Can't send reply: " + reply + ", because there isn't reply topic for original message: " + original);
		}
		replyTopic.sendMessage(reply);
		replyTopicMap.remove(original);
	}

	@Override
	public void replyToMessage(ChipsterMessage original, ChipsterMessage reply, String replyChannel) throws JMSException {
		replyToMessage(original, reply);
	}

	@Override
	public void close() throws JMSException {
		// ignore
	}

	@Override
	public AuthenticationRequestListener getAuthenticationListener() {
		throw new UnsupportedOperationException("not supported by DirectMessagingEndpoint");
	}

	@Override
	public void setAuthenticationListener(
			AuthenticationRequestListener authenticationListener) {
		throw new UnsupportedOperationException("not supported by DirectMessagingEndpoint");
	}

	@Override
	public String getSessionID() {
		throw new UnsupportedOperationException("not supported by DirectMessagingEndpoint");
	}

	@Override
	public void setSessionID(String sessionID) {
		throw new UnsupportedOperationException("not supported by DirectMessagingEndpoint");
	}

	public String getUsername() {
		return username;
	}

	@Override
	public void sendMessageToClientReplyChannel(Destination destination, ChipsterMessage msg)
			throws JMSException {
		throw new UnsupportedOperationException("not supported by DirectMessagingEndpoint");
	}

}
