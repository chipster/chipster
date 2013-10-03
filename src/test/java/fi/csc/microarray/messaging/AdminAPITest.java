package fi.csc.microarray.messaging;

import javax.jms.JMSException;

import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

import fi.csc.microarray.messaging.MessagingTopic.AccessMode;

public class AdminAPITest extends MessagingTestBase {

	@Before
	public void setUp() throws Exception {
		super.setUp();
	}

	@Test
	public void testAdminAPI() throws JMSException, InterruptedException {
		
		MessagingTopic adminTopic = this.endpoint.createTopic(Topics.Name.ADMIN_TOPIC, AccessMode.READ_WRITE);
		AdminAPI api = new AdminAPI(adminTopic, null);
		Assert.assertTrue(api.areAllServicesUp(true));
		
	}
}
