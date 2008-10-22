/*
 * Created on Jan 28, 2005
 *
 */
package fi.csc.microarray.messaging;

/**
 * Base class for implementing MonitoredNode. This class is thread safe.
 * Implemented by calling requestProcessingStarted and requestProcessingStopped
 * appropriately.
 * 
 * @see #requestProcessingStarted()
 * @see #requestProcessingStopped(long)
 * @author akallio
 */
public abstract class MonitoredNodeBase extends NodeBase implements MonitoredNode {

	volatile private long requestsInProcessing = 0;
	volatile private long lastProcessingTime = 0;
	
	public long countRequestsInProcessing() {
		return requestsInProcessing;
	}

	public long getLastProcessingTime() {
		return lastProcessingTime;		
	}

	/**
	 * Subclass must call this when new request is received.
	 * 
	 * @return start time, caller must save this
	 */
	protected long requestProcessingStarted() {
		requestsInProcessing++;
		return System.currentTimeMillis();
	}

	/**
	 * Subclass must call this when new request is processed.
	 * @param startTime time returned by requestProcessingStarted()
	 * @see fi.csc.microarray.messaging.MonitoredNodeBase#requestProcessingStarted()
	 */
	protected void requestProcessingStopped(long startTime) {
		requestsInProcessing--;
		lastProcessingTime = System.currentTimeMillis() - startTime;		
	}
}
