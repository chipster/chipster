package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.util.Queue;
import java.util.concurrent.BlockingQueue;

import javax.swing.SwingUtilities;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;


public abstract class SingleThreadAreaRequestHandler extends AreaRequestHandler {

	private BlockingQueue<AreaRequest> areaRequestQueue;

	private boolean poison = false;

	public SingleThreadAreaRequestHandler(Queue<AreaRequest> areaRequestQueue, AreaResultListener areaResultListener) {

		super(areaRequestQueue, areaResultListener);
		this.areaRequestQueue = (BlockingQueue<AreaRequest>) areaRequestQueue;
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
						if ((areaRequest = areaRequestQueue.take()) != null) {
							areaRequest.getStatus().areaRequestCount = areaRequestQueue.size();
							processAreaRequest(areaRequest);
						}
					} catch (InterruptedException e) {
						e.printStackTrace();
					}
				}
			}
		};

		thread.setDaemon(true);
		thread.setName(getClass().getSimpleName());
		thread.start();
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
	
	@Override
	public void setQueue(Queue<AreaRequest> queue) {
		this.areaRequestQueue = (BlockingQueue<AreaRequest>) queue;
	}
}
