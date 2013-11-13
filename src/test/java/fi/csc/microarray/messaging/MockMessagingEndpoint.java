package fi.csc.microarray.messaging;

import java.util.HashMap;

import javax.jms.JMSException;

import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.Topics.Name;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.messaging.message.ChipsterMessage;

public class MockMessagingEndpoint implements MessagingEndpoint {

	private HashMap<ChipsterMessage, MockMessagingTopic> replyTopicMap = new HashMap<ChipsterMessage, MockMessagingTopic>();
	private HashMap<Name, MockMessagingTopic> topicMap = new HashMap<Name, MockMessagingTopic>();
	
	public void mapReplyTopic(ChipsterMessage original, MockMessagingTopic replyTopic) {
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
			topicMap.put(topicName, new MockMessagingTopic(topicName.name(), accessMode, this));
		}
		
		return topicMap.get(topicName);
	}

	@Override
	public void replyToMessage(ChipsterMessage original, ChipsterMessage reply)	throws JMSException {		
		replyTopicMap.get(original).sendMessage(reply);
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
		return null; // not implemented
	}

	@Override
	public void setAuthenticationListener(
			AuthenticationRequestListener authenticationListener) {
		// ignore
	}

	@Override
	public String getSessionID() {
		return null; // not implemented
	}

	@Override
	public void setSessionID(String sessionID) {
		// ignore
	}

}
