package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.util.Queue;
import java.util.concurrent.LinkedBlockingDeque;

import javax.swing.SwingUtilities;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;


public abstract class SingleThreadAreaRequestHandler extends AreaRequestHandler {

	private LinkedBlockingDeque<AreaRequest> areaRequestQueue;

	private boolean poison = false;

	private Region dataRegion;

	public SingleThreadAreaRequestHandler(Queue<AreaRequest> areaRequestQueue, AreaResultListener areaResultListener) {

		super(areaRequestQueue, areaResultListener);
		this.areaRequestQueue = (LinkedBlockingDeque<AreaRequest>) areaRequestQueue;
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
					try {
						
						if ((areaRequest = areaRequestQueue.takeFirst()) != null) {
							areaRequest.getStatus().areaRequestCount = areaRequestQueue.size();
							
							synchronized (this) {								
								if (dataRegion == null || (dataRegion != null && dataRegion.intersects(areaRequest))) {
									processAreaRequest(areaRequest);
								} else {
									//skip this request, because the data isn't needed anymore
								}							
							}
						}
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
				clean();
			}

		};

		thread.setDaemon(true);
		thread.setName(getClass().getSimpleName());
		thread.start();
	}

	/**
	 * Override this method to close files etc.
	 */
	public void clean() {		
	}

	protected void processAreaRequest(AreaRequest areaRequest) {

		if (areaRequest.getStatus().poison) {

			this.areaResultListener = null;
			this.poison = true;
			
			clean();
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
	
	@Override
	public void setQueue(Queue<AreaRequest> queue) {
		this.areaRequestQueue = (LinkedBlockingDeque<AreaRequest>) queue;
	}

	public boolean hasNewRequest() {
		return areaRequestQueue.size() > 0;
	}

	public void setDataRegion(Region dataRegion) {
		synchronized (this) {			
			this.dataRegion = dataRegion;
		}
	}

	public Region getDataRegion() {
		return dataRegion;
	}
}
