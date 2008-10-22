/*
 * Created on Mar 1, 2005
 *
 */
package fi.csc.microarray.client.tasks;

/**
 * @author Aleksi Kallio
 *
 */
public class TaskException extends Exception {
	
    public TaskException(String msg) {
        super(msg);
    }
    
	public TaskException(Exception e) {
		this(e.getMessage(), e);
	}
	
	public TaskException(String msg, Throwable t) {
		super(msg, t);
	}

}
