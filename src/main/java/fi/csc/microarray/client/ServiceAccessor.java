package fi.csc.microarray.client;

import java.util.Collection;

import fi.csc.microarray.client.operation.OperationCategory;
import fi.csc.microarray.client.tasks.TaskExecutor;
import fi.csc.microarray.databeans.DataManager;
import fi.csc.microarray.messaging.auth.AuthenticationRequestListener;

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

	public String checkRemoveServices() throws Exception;

	public void fetchDescriptions(String module) throws Exception;

	public Collection<OperationCategory> getVisibleCategories() throws Exception;

	public Collection<OperationCategory> getHiddenCategories() throws Exception;

	public TaskExecutor getTaskExecutor() throws Exception;

}
