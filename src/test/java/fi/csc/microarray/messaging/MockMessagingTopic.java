package fi.csc.microarray.messaging;

import java.util.LinkedList;

import javax.jms.JMSException;

public class MockMessagingTopic extends MessagingTopic {

	private LinkedList<MessagingListener> listeners = new LinkedList<MessagingListener>();
	private MockMessagingEndpoint mockEndpoint;
	
	public MockMessagingTopic(String topicName, AccessMode accessMode, MockMessagingEndpoint endpoint)
			throws JMSException {
		super(null, topicName, Type.MOCK, AccessMode.WRITE, endpoint);
		this.mockEndpoint = endpoint;
	}
	
	@Override
	public void sendMessage(fi.csc.microarray.messaging.message.ChipsterMessage message) throws JMSException {
		// MockMessagingEndpoint guarantees that everyone has the same MockMessagingTopic instance
		for (MessagingListener listener : listeners) {
			listener.onChipsterMessage(message);
		}
	}
	
	@Override
	protected void sendReplyableMessage(fi.csc.microarray.messaging.message.ChipsterMessage message, TempTopicMessagingListener replyListener, MessagingListener authenticationListener) throws JMSException {
		sendReplyableMessage(message, replyListener); // ignore authentication
	}
		
	@Override	
	public void setListener(MessagingListener listener) throws JMSException {
		this.listeners.add(listener);
	}
	
	@Override	
	public void sendReplyableMessage(fi.csc.microarray.messaging.message.ChipsterMessage message, TempTopicMessagingListener replyListener) throws JMSException {
		MockMessagingTopic replyTopic = new MockMessagingTopic("temp topic", AccessMode.READ_WRITE, mockEndpoint);
		replyTopic.setListener(replyListener);
		this.mockEndpoint.mapReplyTopic(message, replyTopic);
		sendMessage(message);
		
	}
	
	@Override		
	public void removeListener() throws JMSException {
		// ignore
	}


}
