package fi.csc.microarray.client.visualisation.methods.gbrowser.message;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;
import java.util.concurrent.LinkedBlockingQueue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.SingleThreadAreaRequestHandler;

/**
 * Collects and resends area results. Used by the {@link GBrowserView} objects to manage incoming area results.
 * 
 * @author Petri Klemelä
 *
 */
public class QueueManager implements AreaResultListener {
	
	private class QueueContext {
		public Queue<AreaRequest> queue;
		public Collection<AreaResultListener> listeners = new ArrayList<AreaResultListener>();
		public AreaRequestHandler requestHandler;
	}

	private Map<AreaRequestHandler, QueueContext> queues = new HashMap<AreaRequestHandler, QueueContext>();

	private QueueContext createQueue(AreaRequestHandler areaRequestHandler) {

		if (!queues.containsKey(areaRequestHandler)) {
			QueueContext context = new QueueContext();
			
			if (areaRequestHandler instanceof SingleThreadAreaRequestHandler) {
				
				context.queue = new LinkedBlockingQueue<AreaRequest>();				
			} else {
				context.queue = new ConcurrentLinkedQueue<AreaRequest>();
			}
			
			try {
			    // create a thread which is an instance of class which is passed
			    // as data fetcher to this method
				areaRequestHandler.setQueue(context.queue);
				areaRequestHandler.setAreaResultListener(this);
				
				context.requestHandler = areaRequestHandler;
				queues.put(areaRequestHandler, context);
				
				if (context.requestHandler.isAlive()) {
					System.err.println(
							"Thread '" + context.requestHandler + "' is poisoned, but still alive. " +
									"A new thread will be started for the upcoming requests.");
					context.requestHandler = (AreaRequestHandler) context.requestHandler.clone();
				} 
				
				context.requestHandler.runThread();
				
				return context;

			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return null;
	}
	
	/**
	 * Remove queue for the given data source.
	 * 
	 * @param file
	 */
	public void removeQueue(DataSource file) {
	    queues.remove(file);
	}
	
	public void addAreaRequest(AreaRequestHandler areaRequestHandler, AreaRequest req, boolean clearQueues) {
		
		req.getStatus().areaRequestHandler = areaRequestHandler;
		QueueContext context = queues.get(areaRequestHandler);

		req.getStatus().maybeClearQueue(context.queue);
		context.queue.add(req);
		
		if (context.requestHandler != null) {
		
			//BlockinQueue takes of care of notifyin of SingleThreadAreaRequestHandlers
			if (!(context.requestHandler instanceof SingleThreadAreaRequestHandler)) {
				context.requestHandler.notifyAreaRequestHandler();
			}
		}		
	}

	public void addResultListener(AreaRequestHandler areaRequestHandler, AreaResultListener listener) {
		
		QueueContext qContext = queues.get(areaRequestHandler);
		if (qContext == null) {
			qContext = createQueue(areaRequestHandler);
		}
		qContext.listeners.add(listener);
	}

	public void processAreaResult(AreaResult areaResult) {

		for (AreaResultListener listener : queues.get(areaResult.getStatus().areaRequestHandler).listeners) {
			listener.processAreaResult(areaResult);
		}
	}

	public void poisonAll() {
		
		for (Entry<AreaRequestHandler, QueueContext> entry : queues.entrySet()) {
			
			DataRetrievalStatus status = new DataRetrievalStatus();
			status.poison = true;
			AreaRequest request = new AreaRequest(new Region(), null, status);
						
			QueueContext context = entry.getValue();
			context.queue.add(request);
			context.requestHandler.notifyAreaRequestHandler();		
		}
	}
}
