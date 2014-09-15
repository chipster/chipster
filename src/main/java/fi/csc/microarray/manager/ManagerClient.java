package fi.csc.microarray.manager;

import java.net.URL;
import java.util.List;

import javax.jms.JMSException;

import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;

// empty impl.
public class ManagerClient {
	
	public ManagerClient(MessagingEndpoint endpoint) throws JMSException {
		if (endpoint != null) {
			endpoint.createTopic(Topics.Name.JOB_LOG_TOPIC, AccessMode.WRITE);
		}
	}
	
	public void urlRequest(String username, URL url) {
		// url should be secured in some way before transmitting
		
	}

	@Deprecated
	public void publicUrlRequest(String username, URL url) {
		
	}

	public void publicFilesRequest(String username, List<URL> files) {
		
	}
}
