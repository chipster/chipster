package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.FsfStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

/**
 * Collects and resends area results. Used by the {@link View} objects to manage incoming area results.
 * 
 * @author Petri Klemelä
 *
 */
public class QueueManager implements AreaResultListener {

	private class QueueContext {
		public Queue<AreaRequest> queue;
		public Collection<AreaResultListener> listeners = new ArrayList<AreaResultListener>();
		public AreaRequestHandler thread;
	}

	private Map<DataSource, QueueContext> queues = new HashMap<DataSource, QueueContext>();

	private QueueContext createQueue(DataSource file) {

		if (!queues.containsKey(file)) {
			QueueContext context = new QueueContext();
			context.queue = new ConcurrentLinkedQueue<AreaRequest>();
			try {
			    // create a thread which is an instance of class which is passed
			    // as data fetcher to this method
				context.thread = file.getRequestHandler().getConstructor(DataSource.class,
				        Queue.class, AreaResultListener.class).
				        newInstance(file, context.queue, this);

				queues.put(file, context);
				context.thread.start();
				
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

	public void addAreaRequest(DataSource file, AreaRequest req, boolean clearQueues) {
		req.status.file = file;
		QueueContext context = queues.get(file);

		req.status.maybeClearQueue(context.queue);
		context.queue.add(req);
		
		if (context.thread != null) {
			context.thread.notifyAreaRequestHandler();
		}
	}

	public void addResultListener(DataSource file, AreaResultListener listener) {
		
		QueueContext qContext = queues.get(file);
		if (qContext == null) {
			qContext = createQueue(file);
		}
		qContext.listeners.add(listener);
	}

	public void processAreaResult(AreaResult areaResult) {

		for (AreaResultListener listener : queues.get(areaResult.getStatus().file).listeners) {
			listener.processAreaResult(areaResult);
		}
	}

	public void poisonAll() {
		
		for (Entry<DataSource, QueueContext> entry : queues.entrySet()) {
			
			FsfStatus status = new FsfStatus();
			status.poison = true;
			AreaRequest request = new AreaRequest(new Region(), null, status);
						
			QueueContext context = entry.getValue();
			context.queue.add(request);
			context.thread.notifyAreaRequestHandler();		
		}
	}
}
