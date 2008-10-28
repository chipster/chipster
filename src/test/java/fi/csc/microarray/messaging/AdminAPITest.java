package fi.csc.microarray.messaging;

import javax.jms.JMSException;

import org.testng.Assert;
import org.testng.annotations.BeforeSuite;
import org.testng.annotations.Test;

import fi.csc.microarray.AdminAPI;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;

public class AdminAPITest extends MessagingTestBase {

	@BeforeSuite
	public void setUp() throws Exception {
		super.setUp();
	}

	@Test(groups = { "smoke"})
	public void testAdminAPI() throws JMSException, InterruptedException {
		
		MessagingTopic adminTopic = this.endpoint.createTopic(Topics.Name.ADMIN_TOPIC, AccessMode.READ_WRITE);
		AdminAPI api = new AdminAPI(adminTopic, null);
		Assert.assertTrue(api.areAllServicesUp(true));
		
	}
}
