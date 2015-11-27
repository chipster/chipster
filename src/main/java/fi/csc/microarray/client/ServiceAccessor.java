package fi.csc.microarray.client;

import java.util.Collection;

import fi.csc.microarray.client.operation.ToolModule;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.filebroker.FileBrokerClient;
import fi.csc.microarray.messaging.SourceMessageListener;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;
import fi.csc.microarray.messaging.message.FeedbackMessage;
import fi.csc.microarray.module.Module;

/**
 * Interface that client uses to access services of the underlying infrastructure.
 * Abstracts that fact that although services typically are remote, they can also
 * be local if client is run in standalone mode.
 * 
 * @author Aleksi Kallio
 *
 */
public interface ServiceAccessor {

	public static final String ALL_SERVICES_OK = "ok";
	
	public void initialise(DataManager manager, AuthenticationRequestListener authenticationRequestListener) throws Exception;

	public String checkRemoteServices() throws Exception;

	/**
	 * 
	 * @param primaryModule
	 * @return error messages if fetching involved errors, otherwise empty string
	 * @throws Exception
	 */
	public String fetchDescriptions(Module primaryModule) throws Exception;

	public TaskExecutor getTaskExecutor();

	public void close();

	public SourceMessageListener retrieveSourceCode(String id) throws Exception;

	public FileBrokerClient getFileBrokerClient();

	public void sendFeedbackMessage(FeedbackMessage message) throws Exception;

	public Collection<ToolModule> getModules();
	
	public boolean isStandalone();

}
