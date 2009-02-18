package fi.csc.microarray.manager;

import java.net.URL;

import javax.jms.JMSException;

import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;

// empty impl.
public class ManagerClient {
	
	public ManagerClient(MessagingEndpoint endpoint) throws JMSException {
		endpoint.createTopic(Topics.Name.MANAGER_TOPIC, AccessMode.WRITE);
	}
	
	public void urlRequest(String username, URL url) {
		// url should be secured in some way before transmitting
		
	}

}
