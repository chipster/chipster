package fi.csc.microarray.messaging;

import javax.jms.JMSException;

import org.junit.Assert;
import org.junit.Test;

import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.Topics.Name;
import fi.csc.microarray.messaging.message.ChipsterMessage;

public class DirectMessagingEndpointTest {
	
	private static final String USERNAME = "authenticated-username";
	
	private ChipsterMessage receivedMessageTopic1;
	private ChipsterMessage receivedMessageTopic2;
	private ChipsterMessage receivedReplyMessage;

	@Test
	public void test() throws JMSException {
		DirectMessagingEndpoint endpoint = new DirectMessagingEndpoint();
		DirectMessagingEndpoint authenticatingEndpoint = new DirectMessagingEndpoint(USERNAME);
		
		testEndpoint(endpoint, null);
		testEndpoint(authenticatingEndpoint, USERNAME);
	}
		
	public void testEndpoint(MessagingEndpoint endpoint, String username) throws JMSException {
		
		MessagingTopic topic1 = endpoint.createTopic(Name.TEST_TOPIC, AccessMode.WRITE);
		MessagingTopic topic2 = endpoint.createTopic(Name.REQUEST_TOPIC, AccessMode.WRITE);
		
		topic1.setListener(new MessagingListener() {

			@Override
			public void onChipsterMessage(ChipsterMessage msg) {
				receivedMessageTopic1 = msg;
			}			
		});
		
		topic2.setListener(new MessagingListener() {

			@Override
			public void onChipsterMessage(ChipsterMessage msg) {
				receivedMessageTopic2 = msg;
			}			
		});
		
		ChipsterMessage message = new ChipsterMessage() { };
		ChipsterMessage replyMessage = new ChipsterMessage() { };
		
		//send a simple message
		
		topic1.sendMessage(message);		
		Assert.assertEquals(message.getMessageID(), receivedMessageTopic1.getMessageID());
		Assert.assertNull(receivedMessageTopic2);
		Assert.assertEquals(username, receivedMessageTopic1.getUsername());
		
		//reply to non-replyable message
		
		ChipsterMessage replyTo = receivedMessageTopic1;
		receivedMessageTopic1 = null;
		
		try {
			endpoint.replyToMessage(replyTo, replyMessage);
			Assert.fail();
		} catch (Exception e) {			
		}
		
		Assert.assertNull(receivedMessageTopic1);
		Assert.assertNull(receivedMessageTopic2);
		Assert.assertNull(receivedReplyMessage);
		
		//send replyable message
		
		topic1.sendReplyableMessage(message, new TempTopicMessagingListener() {
			
			@Override
			public void onChipsterMessage(ChipsterMessage msg) {
				receivedReplyMessage = msg;
			}
			
			@Override
			public void setTempTopic(MessagingTopic tempTopic) {
			}
			
			@Override
			public void cleanUp() {
			}

			@Override
			public void cancel() {
			}
		});		
		Assert.assertEquals(message.getMessageID(), receivedMessageTopic1.getMessageID());
		Assert.assertNull(receivedMessageTopic2);
		Assert.assertNull(receivedReplyMessage);
		Assert.assertEquals(username, receivedMessageTopic1.getUsername());
		
		replyTo = receivedMessageTopic1;
		receivedMessageTopic1 = null;
		
		//reply to replyable message
		
		endpoint.replyToMessage(replyTo, replyMessage);
		
		Assert.assertNull(receivedMessageTopic1);
		Assert.assertNull(receivedMessageTopic2);
		Assert.assertEquals(replyMessage.getMessageID(), receivedReplyMessage.getMessageID());
		Assert.assertEquals(username, receivedReplyMessage.getUsername());
		
		receivedReplyMessage = null;
	}
}
