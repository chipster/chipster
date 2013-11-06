package fi.csc.microarray.messaging;

import javax.jms.JMSException;

import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.Topics.Name;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.messaging.message.ChipsterMessage;

public class MockMessagingEndpoint implements MessagingEndpoint {


	@Override
	public MessagingTopic createTopic(Name topicName, AccessMode accessMode)
			throws JMSException {
		return new MockMessagingTopic(topicName.name(), accessMode, this);
	}

	@Override
	public void replyToMessage(ChipsterMessage original, ChipsterMessage reply)
			throws JMSException {
	}

	@Override
	public void replyToMessage(ChipsterMessage original, ChipsterMessage reply,
			String replyChannel) throws JMSException {
	}

	@Override
	public void close() throws JMSException {
	}

	@Override
	public AuthenticationRequestListener getAuthenticationListener() {
		return null;
	}

	@Override
	public void setAuthenticationListener(
			AuthenticationRequestListener authenticationListener) {
	}

	@Override
	public String getSessionID() {
		return null;
	}

	@Override
	public void setSessionID(String sessionID) {
	}

}
