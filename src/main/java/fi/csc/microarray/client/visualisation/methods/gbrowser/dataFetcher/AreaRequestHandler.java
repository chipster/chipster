package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.Queue;

import javax.swing.SwingUtilities;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;

/**
 * The active thread of the processing layer. Receives area requests and sends out
 * file request to the lower data source layer.
 *  
 * @author Petri Klemelä, Aleksi Kallio
 *
 */
public abstract class AreaRequestHandler implements Cloneable {

	private Queue<AreaRequest> areaRequestQueue;
	protected AreaResultListener areaResultListener;

	private boolean poison = false;
	protected Thread thread;

	public AreaRequestHandler(Queue<AreaRequest> areaRequestQueue, AreaResultListener areaResultListener) {

		super();
		this.areaRequestQueue = areaRequestQueue;
		this.areaResultListener = areaResultListener;
	}

	/**
	 * Constantly look for new area requests in the request queue
	 * and pass all incoming requests to request processor.
	 */
	public void runThread() {
		
		thread = new Thread() {
			public synchronized void run() {		

				while (!poison) {
					AreaRequest areaRequest;
					if ((areaRequest = areaRequestQueue.poll()) != null) {
						areaRequest.getStatus().areaRequestCount = areaRequestQueue.size();
						processAreaRequest(areaRequest);
					}

					boolean isWorkToDo = checkOtherQueues();

					if (areaRequest == null && !isWorkToDo) {
						try {
							thread.wait();
						} catch (InterruptedException e) {
							// just poll the queues and wait again
						} catch (IllegalMonitorStateException e) {
							e.printStackTrace();
							poison = true;
						}
					}
				}
			}
		};
		
		thread.setDaemon(true);
		thread.setName(getClass().getSimpleName());
		thread.start();
	}

	/**
	 * @return true if should be called again before putting thread to wait
	 */
	protected boolean checkOtherQueues() {		
		return false; // hook for checking fileFetcherQueue
	}

	public synchronized void notifyAreaRequestHandler() {
		synchronized (thread) {			
			thread.notifyAll();
		}
	}

	protected void processAreaRequest(AreaRequest areaRequest) {

		if (areaRequest.getStatus().poison) {

			this.areaResultListener = null;
			this.poison = true;
		}
	}

	/**
	 * Pass the result to be visualised in GUI.
	 * 
	 * @param areaResult
	 */
	public void createAreaResult(final AreaResult areaResult) {

		SwingUtilities.invokeLater(new Runnable() {
						
			public void run() {
				
				if (areaResultListener != null) {
					areaResultListener.processAreaResult(areaResult);
				}
			}						
		});
	}

	public void setAreaResultListener(AreaResultListener areaResultListener) {
		this.areaResultListener = areaResultListener;
	}

	public void setQueue(Queue<AreaRequest> queue) {
		this.areaRequestQueue = queue;
	}
	
	@Override
	public Object clone() throws CloneNotSupportedException {
		return super.clone();
	}

	public boolean isAlive() {
		if (thread == null) {
			return false;
		} else {
			return thread.isAlive();
		}
	}
}
