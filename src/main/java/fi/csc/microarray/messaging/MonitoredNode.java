package fi.csc.microarray.messaging;

/**
 * Adds possibility to extract state information from the component. This
 * information is used by administrative tools. You can use MonitoredNodeBase 
 * for easier implementation.
 * 
 * @author akallio
 */
// did not extend node, but what was the original reason for that?
public interface MonitoredNode extends Node {
	/**
	 * Returns number of requests currently in processing.
	 */
	public long countRequestsInProcessing();
	
	/**
	 * Returns the processing time in milliseconds for the last successful request.  
	 */
	public long getLastProcessingTime();
	
}
