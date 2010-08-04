package fi.csc.microarray.client;

import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingTopic;

public class RemoteServiceAccessor implements ServiceAccessor {

	protected TaskExecutor taskExecutor;
	protected MessagingEndpoint endpoint;
	protected MessagingTopic requestTopic;

}
