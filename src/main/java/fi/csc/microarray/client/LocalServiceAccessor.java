package fi.csc.microarray.client;

import java.util.Collection;
import java.util.LinkedList;

import javax.jms.JMSException;

import fi.csc.microarray.client.operation.OperationCategory;
import fi.csc.microarray.client.tasks.LocalTaskExecutor;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.messaging.SourceMessageListener;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.messaging.message.FeedbackMessage;

public class LocalServiceAccessor implements ServiceAccessor {

	private DataManager manager;

	@Override
	public String checkRemoveServices() throws Exception {
		return ServiceAccessor.ALL_SERVICES_OK;
	}

	@Override
	public void close() throws Exception {
		// do nothing
	}

	@Override
	public void fetchDescriptions(String module) throws Exception {
		// do nothing
	}

	@Override
	public FileBrokerClient getFileBrokerClient() throws Exception {
		throw new UnsupportedOperationException("not supported in standalone mode");
	}

	@Override
	public Collection<OperationCategory> getHiddenCategories() {
		return new LinkedList<OperationCategory>();
	}

	@Override
	public TaskExecutor getTaskExecutor() {
		try {
			return new LocalTaskExecutor(manager);
		} catch (JMSException e) {
			// TODO Auto-generated catch block
			throw new RuntimeException(e);
		}
	}

	@Override
	public Collection<OperationCategory> getVisibleCategories() {
		return new LinkedList<OperationCategory>();
	}

	@Override
	public void initialise(DataManager manager, AuthenticationRequestListener authenticationRequestListener) throws Exception {
		this.manager = manager;
		// we are not interested in authenticationRequestListener
	}

	@Override
	public SourceMessageListener retrieveSourceCode(String id) throws Exception {
		throw new UnsupportedOperationException("not supported in standalone mode");
	}

	@Override
	public void sendFeedbackMessage(FeedbackMessage message) throws Exception {
		throw new UnsupportedOperationException("not supported in standalone mode");
	}

}
