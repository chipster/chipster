/*
 * Created on Feb 10, 2005
 *
 */
package fi.csc.microarray.client.tasks;

/**
 * Listener for events related to a single task. Guaranteed to be called in Event Dispatch Thread.
 * 
 * @author Aleksi Kallio
 *
 */
public interface TaskEventListener {
	
	public void onStateChange(Task task, Task.State oldState, Task.State newState);
}
