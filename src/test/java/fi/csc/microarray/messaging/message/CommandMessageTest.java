package fi.csc.microarray.messaging.message;

import javax.jms.JMSException;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import fi.csc.microarray.messaging.MessagingListener;
import fi.csc.microarray.messaging.MessagingTestBase;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.Topics;

public class CommandMessageTest extends MessagingTestBase {

	
	private String command;
	
	@Before
	public void setUp() throws Exception {
		super.setUp();
	}
	
	@Test
	public void testSend() throws JMSException {
		CommandMessage msg = new CommandMessage("test");
		endpoint.createTopic(Topics.Name.TEST_TOPIC, AccessMode.WRITE).sendMessage(msg);
	}

	@Test
	public void testReceive() throws JMSException, InterruptedException {
		
		endpoint.createTopic(Topics.Name.TEST_TOPIC, AccessMode.READ).setListener(new MessagingListener() {
			public void onChipsterMessage(ChipsterMessage msg) {
				CommandMessage cmsg = (CommandMessage)msg;
				command = cmsg.getCommand();
			}
		});
		testSend();
		Thread.sleep(1000); // is this a long enough pause?
		Assert.assertTrue(command.equals("test"));
	}
}
