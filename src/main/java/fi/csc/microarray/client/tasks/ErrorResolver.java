package fi.csc.microarray.client.tasks;

/**
 * 
 * @author Aleksi Kallio
 *
 */
public interface ErrorResolver {
	
	public Task.State resolveErrorState(Task task);
	public String resolveErrorMessage(Task task);
}
