package fi.csc.microarray.client.visualisation.methods.gbrowser.stack;

import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.LinkedBlockingDeque;

import javax.swing.SwingUtilities;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRetrievalStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;


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
						
						//report before waiting in queue.take()
						reportQueueSize(false);
						
						if ((areaRequest = areaRequestQueue.takeFirst()) != null) {
							
							//report after wait
							reportQueueSize(true);
							
							//areaRequest.getStatus().areaRequestCount = areaRequestQueue.size();
							
							synchronized (this) {
								

								if (dataRegion == null || 
										areaRequest.getStatus().poison || //poison requests would cause nullPointer on the next line
										areaRequest instanceof GeneRequest || //searched gene may be in other chromosome
										(dataRegion != null && dataRegion.intersects(areaRequest))) {
									
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

	/**
	 * This background thread is processing area requests in the order they appear from the queue. 
	 * Sometimes the queue is so long that some request aren't needed anymore. This method
	 * bypasses the queue and sets the current view region, so that old requests are removed, if those don't intercept
	 * with this region.
	 * 
	 * @param dataRegion
	 */
	public void setDataRegion(Region dataRegion) {
		synchronized (this) {			
			if (dataRegion != null && dataRegion.start != null && dataRegion.end != null) {
				//Create a new Region instance for this thread
				Region region = new Region((long)dataRegion.start.bp, (long)dataRegion.end.bp, new Chromosome(dataRegion.start.chr));
				this.dataRegion = region;
			}
		}
	}

	public Region getDataRegion() {
		return dataRegion;
	}

	/**
	 * @param increaseByOne set true to indicate that queue size should be increased by one to count the
	 * request which is taken from the queue, but not processed yet
	 */
	private void reportQueueSize(boolean increaseByOne) {
		
		
		DataRetrievalStatus status = new DataRetrievalStatus();
		status.areaRequestHandler = this;
		status.areaRequestCount = areaRequestQueue.size() + (increaseByOne ? 1 : 0);		
		List<RegionContent> emptyList = new LinkedList<RegionContent>();
		createAreaResult(new AreaResult(status, emptyList));		
	}
}
