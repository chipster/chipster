package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Map.Entry;
import java.util.Queue;
import java.util.concurrent.LinkedBlockingDeque;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResultListener;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataStatus;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;

/**
 * Collects and resends dataResults. Used by the {@link GBrowserView} objects to manage incoming dataResults.
 * 
 * @author Petri Klemelä
 *
 */
public class QueueManager {
	
	private class QueueContext {
		public Queue<DataRequest> queue;
		public Collection<DataResultListener> listeners = new ArrayList<DataResultListener>();
		public DataThread dataThread;
	}

	private Map<DataThread, QueueContext> queues = new HashMap<DataThread, QueueContext>();
	private GBrowserView view;

	public QueueManager(GBrowserView view) {
		this.view = view;
	}

	private QueueContext createQueue(DataThread dataThread) {

		if (!queues.containsKey(dataThread)) {
			QueueContext context = new QueueContext();
			
			context.queue = new LinkedBlockingDeque<DataRequest>();				
			
			try {

				dataThread.setQueue(context.queue);
				dataThread.setQueueManager(this);
				
				context.dataThread = dataThread;
				queues.put(dataThread, context);
				
				if (context.dataThread.isAlive()) {
					System.err.println(
							"Thread '" + context.dataThread + "' is poisoned, but still alive. " +
									"A new thread will be started for the upcoming requests.");
					context.dataThread = (DataThread) context.dataThread.clone();
				} 
				
				context.dataThread.runThread();
				
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
	
	public void addDataRequest(DataThread dataThread, DataRequest req, Region dataRegion) {
		
		req.getStatus().setDataThread(dataThread);
		QueueContext context = queues.get(dataThread);

		context.dataThread.setDataRegion(dataRegion);
		
		context.queue.add(req);		
	}

	public void addDataResultListener(DataThread dataThread, DataResultListener listener) {
		
		QueueContext qContext = queues.get(dataThread);
		if (qContext == null) {
			qContext = createQueue(dataThread);
		}
		qContext.listeners.add(listener);
	}

	public void processDataResult(DataResult dataResult) {
		
		for (DataResultListener listener : queues.get(dataResult.getStatus().getDataThread()).listeners) {
			
//			long t = System.currentTimeMillis();
			
			if (listener instanceof Track) {
				Track track = (Track) listener;
				
				if (track.isSuitableViewLength()) {
					listener.processDataResult(dataResult);	
				}
				
			} else {
				listener.processDataResult(dataResult);
			}
						
//			t = System.currentTimeMillis() - t;
//			
//			if (t > 1) {
//				System.out.println(listener + "\t" + t);
//			}					
		}
		
		view.redraw();
	}

	public void poisonAll() {
		
		for (Entry<DataThread, QueueContext> entry : queues.entrySet()) {
			
			DataStatus status = new DataStatus();
			status.poison = true;
			DataRequest request = new DataRequest(new Region(), null, status);
						
			QueueContext context = entry.getValue();
			context.queue.add(request);		
		}
	}

	public Collection<DataResultListener> getDataResultListeners(DataResult dataResult) {
		return queues.get(dataResult.getStatus().getDataThread()).listeners;
	}

	public void clearDataResultListeners() {
		for (QueueContext context : queues.values()) {
			context.listeners.clear();
		}
	}
}
