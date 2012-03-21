package fi.csc.microarray.client;

import java.util.Collection;

import javax.jms.JMSException;

import fi.csc.microarray.client.operation.ToolModule;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.filebroker.JMSFileBrokerClient;
import fi.csc.microarray.messaging.AdminAPI;
import fi.csc.microarray.messaging.DescriptionMessageListener;
import fi.csc.microarray.messaging.MessagingEndpoint;
import fi.csc.microarray.messaging.MessagingTopic;
import fi.csc.microarray.messaging.NodeBase;
import fi.csc.microarray.messaging.SourceMessageListener;
import fi.csc.microarray.messaging.Topics;
import fi.csc.microarray.messaging.MessagingTopic.AccessMode;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.messaging.message.CommandMessage;
import fi.csc.microarray.messaging.message.FeedbackMessage;
import fi.csc.microarray.module.Module;

public class RemoteServiceAccessor implements ServiceAccessor {

	protected MessagingEndpoint endpoint;
	protected MessagingTopic requestTopic;
	protected TaskExecutor taskExecutor;
	private FileBrokerClient filebrokerClient;

	private NodeBase nodeSupport = new NodeBase() {
		public String getName() {
			return "client";
		}
	};
	private Collection<ToolModule> modules = null;

	public void initialise(DataManager manager, AuthenticationRequestListener authenticationRequestListener) throws MicroarrayException, JMSException {
		this.endpoint = new MessagingEndpoint(nodeSupport, authenticationRequestListener);
	    this.requestTopic = endpoint.createTopic(Topics.Name.REQUEST_TOPIC,AccessMode.WRITE);
		this.taskExecutor = new TaskExecutor(endpoint, manager);
		this.filebrokerClient = new JMSFileBrokerClient(endpoint.createTopic(Topics.Name.URL_TOPIC, AccessMode.WRITE));
	}

	public TaskExecutor getTaskExecutor() {
		if (taskExecutor == null) {
			throw new IllegalStateException("initialise(...) must be called first");
		}
		return taskExecutor;
	}
	
	@Override
	public String checkRemoteServices() throws Exception {
		AdminAPI api = new AdminAPI(endpoint.createTopic(Topics.Name.ADMIN_TOPIC, AccessMode.READ_WRITE), null);
		if (api.areAllServicesUp(true)) {
			return ALL_SERVICES_OK;
		} else {
			return "required services are not available (" + api.getErrorStatus() + ")";
		}				
	}

	@Override
	public void fetchDescriptions(Module module) throws Exception {
        DescriptionMessageListener descriptionListener = new DescriptionMessageListener(module.getServerModuleNames());
		try {
			this.requestTopic.sendReplyableMessage(new CommandMessage(CommandMessage.COMMAND_DESCRIBE), descriptionListener);
			descriptionListener.waitForResponse();
		} finally {
			descriptionListener.cleanUp();
		}
		this.modules = descriptionListener.getModules();
	}


	@Override
	public Collection<ToolModule> getModules() {
		if (modules == null) {
			throw new IllegalStateException("fetchDescriptions(...) must be called first");
		}
		return modules;
	}
	
	@Override
	public void close() throws Exception {
		endpoint.close();
	}

	@Override
	public SourceMessageListener retrieveSourceCode(String id) throws Exception {
		SourceMessageListener sourceListener = new SourceMessageListener();
		CommandMessage commandMessage = new CommandMessage(CommandMessage.COMMAND_GET_SOURCE);
		commandMessage.addParameter(id);
		this.requestTopic.sendReplyableMessage(commandMessage, sourceListener);
		return sourceListener;
	}

	@Override
	public FileBrokerClient getFileBrokerClient() throws Exception {
        return filebrokerClient;
	}

	@Override
	public void sendFeedbackMessage(FeedbackMessage message) throws Exception {
        MessagingTopic requestTopic = endpoint.createTopic(Topics.Name.FEEDBACK_TOPIC, AccessMode.WRITE);
        requestTopic.sendMessage(message);
	}

	@Override
	public boolean isStandalone() {
		return false;
	}

}
