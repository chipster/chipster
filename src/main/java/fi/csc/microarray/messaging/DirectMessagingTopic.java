package fi.csc.microarray.messaging;

import java.util.LinkedList;

import javax.jms.JMSException;

import fi.csc.microarray.messaging.message.ChipsterMessage;

public class DirectMessagingTopic extends MessagingTopic {

	private LinkedList<MessagingListener> listeners = new LinkedList<MessagingListener>();
	private DirectMessagingEndpoint endpoint;
	
	public DirectMessagingTopic(String topicName, AccessMode accessMode, DirectMessagingEndpoint endpoint)
			throws JMSException {
		super(null, topicName, Type.MOCK, AccessMode.WRITE, endpoint);
		this.endpoint = endpoint;
	}
	
	@Override
	public void sendMessage(fi.csc.microarray.messaging.message.ChipsterMessage message) throws JMSException {
		
		setUsername(message);
		
		// DirectMessagingEndpoint guarantees that everyone has the same DirectMessagingTopic instance
		for (MessagingListener listener : listeners) {
			if (listener instanceof DirectMessagingListener) {
				DirectMessagingListener directListener = (DirectMessagingListener) listener;
				directListener.onChipsterMessage(message, endpoint);
			} else {
				listener.onChipsterMessage(message);
			}
		}
	}
	
	private void setUsername(ChipsterMessage message) {
		if (endpoint.getUsername() != null) {
			message.setUsername(endpoint.getUsername());
		}
	}

	@Override
	protected void sendReplyableMessage(fi.csc.microarray.messaging.message.ChipsterMessage message, TempTopicMessagingListener replyListener, AuthMessagingListener authenticationListener) throws JMSException {
		sendReplyableMessage(message, replyListener); // ignore authentication
	}
		
	@Override	
	public void setListener(MessagingListener listener) throws JMSException {
		this.listeners.add(listener);
	}
	
	@Override	
	public void sendReplyableMessage(fi.csc.microarray.messaging.message.ChipsterMessage message, TempTopicMessagingListener replyListener) throws JMSException {
		DirectMessagingTopic replyTopic = new DirectMessagingTopic("temp topic", AccessMode.READ_WRITE, endpoint);
		replyTopic.setListener(replyListener);
		this.endpoint.mapReplyTopic(message, replyTopic);
		setUsername(message);
		sendMessage(message);
		
	}
	
	@Override		
	public void removeListener() throws JMSException {
		throw new UnsupportedOperationException("not supported by DirectMessagingTopic");
	}
}
