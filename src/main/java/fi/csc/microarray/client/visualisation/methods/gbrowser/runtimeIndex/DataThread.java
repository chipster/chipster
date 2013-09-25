package fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex;

import java.util.LinkedList;
import java.util.List;
import java.util.Queue;
import java.util.concurrent.LinkedBlockingDeque;

import javax.swing.SwingUtilities;

import fi.csc.microarray.client.visualisation.methods.gbrowser.GBrowser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.QueueManager;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.GBrowserException;


/**
 * This class contains a generic implementations for managing the execution and messaging of the background
 * data threads. These threads extend this class and usually convert String or binary content of the file to Java objects.  Therefore names 
 * of those classes usually include a word 'Conversion'. 
 * 
 * @author klemela
 */
public abstract class DataThread {

	private LinkedBlockingDeque<DataRequest> dataRequestQueue;
	protected QueueManager queueManager;
	protected Thread thread;


	private boolean poison = false;

	private Region dataRegion;

	private GBrowser browser;
	private DataSource dataSource;

	public DataThread(GBrowser browser, DataSource dataSource) {
		
		this.browser = browser;
		this.dataSource = dataSource;
	}

	/**
	 * Constantly look for new data requests in the request queue
	 * and pass all incoming requests to request processor.
	 */
	public void runThread() {

		thread = new Thread() {
			public synchronized void run() {		

				while (!poison) {
					DataRequest dataRequest;
					try {
						
						//report before waiting in queue.take()
						reportQueueSize(false);
						
						if ((dataRequest = dataRequestQueue.takeFirst()) != null) {
							
							//report queue length after wait
							reportQueueSize(true);						
							
							boolean processRequest;
							
							synchronized (this) {
								
								processRequest = 
										dataRegion == null || 
										dataRequest.getStatus().poison || //poison requests would cause nullPointer on the next line
										dataRequest instanceof GeneRequest || //searched gene may be in other chromosome
										(dataRegion != null && dataRegion.intersects(dataRequest));
							}
							
							if (processRequest) {
								
								if (!checkPoison(dataRequest)) {
									try {
										processDataRequest(dataRequest);
									} catch (GBrowserException e) {
										reportException(e);
										poison = true;
									}
								}
							} else {
								//skip this request, because the data isn't needed anymore
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
	
	private boolean checkPoison(DataRequest dataRequest) {
		if (dataRequest.getStatus().poison) {

			this.queueManager = null;
			this.poison = true;
			
			clean();
		}
		
		return dataRequest.getStatus().poison;
	}

	protected abstract void processDataRequest(DataRequest dataRequest) throws GBrowserException;
	/**
	 * Pass the result to be visualised in GUI.
	 * 
	 * @param dataResult
	 */
	public void createDataResult(final DataResult dataResult) {

		SwingUtilities.invokeLater(new Runnable() {
						
			public void run() {
				
				if (queueManager != null) {
					queueManager.processDataResult(dataResult);
				}
			}						
		});
	}	
	
	public void setQueue(Queue<DataRequest> queue) {
		this.dataRequestQueue = (LinkedBlockingDeque<DataRequest>) queue;
	}

	public boolean hasNewRequest() {
		return dataRequestQueue.size() > 0;
	}

	/**
	 * This background thread is processing data requests in the order they appear from the queue. 
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
		
		
		DataStatus status = new DataStatus();
		status.setDataThread(this);
		status.setDataRequestCount(dataRequestQueue.size() + (increaseByOne ? 1 : 0));		
		List<Feature> emptyList = new LinkedList<Feature>();		
		createDataResult(new DataResult(status, emptyList));		
	}

	public void setQueueManager(QueueManager queueManager) {
		this.queueManager = queueManager;
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
	
	private void reportException(final Exception e) {
		SwingUtilities.invokeLater(new Runnable() {
			
			@Override
			public void run() {
				browser.reportException(e);
			}
		});
	}

	public void setDataSource(DataSource dataSource) {
		this.dataSource = dataSource;
	}

	public DataSource getDataSource() {
		return dataSource;
	}
}
