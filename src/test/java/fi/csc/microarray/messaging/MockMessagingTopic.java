package fi.csc.microarray.messaging;

import javax.jms.JMSException;

public class MockMessagingTopic extends MessagingTopic {

	public MockMessagingTopic(String topicName, AccessMode accessMode, MessagingEndpoint endpoint)
			throws JMSException {
		super(null, topicName, Type.MOCK, AccessMode.WRITE, endpoint);
	}
	
	@Override
	public void sendMessage(fi.csc.microarray.messaging.message.ChipsterMessage message) throws JMSException {
		// empty
	}
	
	@Override
	protected void sendReplyableMessage(fi.csc.microarray.messaging.message.ChipsterMessage message, TempTopicMessagingListener replyListener, MessagingListener authenticationListener) throws JMSException {
		// empty
	}
		
	@Override	
	public void setListener(MessagingListener listener) throws JMSException {
		// empty
	}
	
	@Override	
	public void sendReplyableMessage(fi.csc.microarray.messaging.message.ChipsterMessage message, TempTopicMessagingListener replyListener) throws JMSException {
		// empty		
	}
	
	@Override		
	public void removeListener() throws JMSException {
		// empty
	}


}
