/*
 * Created on Jan 28, 2005
 *
 */
package fi.csc.microarray.messaging;

import javax.jms.ExceptionListener;

/**
 * Interface for callback from MessagingEndpoint. Use NodeBase for easier implementation. 
 * 
 * @author akallio
 */
public interface Node extends ExceptionListener {

	/**
	 * Returns name of the node. Used only for informative purposes, so any non-null String is valid.
	 */
	public String getName();
	
	/**
	 * Returns the hostname that the component is run on. 
	 */
	public String getHost();
}
