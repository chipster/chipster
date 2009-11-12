package fi.csc.microarray.client.visualisation.methods.genomeBrowser.dataFetcher;

import java.io.File;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.Map;
import java.util.Queue;
import java.util.concurrent.ConcurrentLinkedQueue;

import fi.csc.microarray.client.visualisation.methods.genomeBrowser.fileFormat.ReadInstructions;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.genomeBrowser.message.AreaResult;

public class QueueManager implements AreaResultListener {

	private class QueueContext {
		public Queue<AreaRequest> queue;
		public Collection<AreaResultListener> listeners = new ArrayList<AreaResultListener>();
		public AreaRequestHandler thread;		
	}

	Map<File, QueueContext> queues = new HashMap<File, QueueContext>(); 

	public void createQueue(File file, ReadInstructions<?> readInstructions){
		createQueue(file, StraightforwardFileParser.class, readInstructions);
	}		

	public void createQueue(File file, Class<? extends AreaRequestHandler> dataFetcher, ReadInstructions<?> readInstructions){
		if(!queues.containsKey(file)){
			QueueContext context = new QueueContext();
			context.queue = new ConcurrentLinkedQueue<AreaRequest>();		
			try {
				context.thread = dataFetcher.getConstructor(
						File.class, 
						Queue.class, 
						AreaResultListener.class, 
						ReadInstructions.class).
						
					newInstance(file, context.queue, this, readInstructions);

				queues.put(file, context);
				context.thread.start();
				
			} catch (Exception e){
				e.printStackTrace();
			}		
		}
	}

	public void addAreaRequest(File file, AreaRequest req, boolean clearQueues){
		req.status.file = file;
		QueueContext context = queues.get(file);
		
		//System.out.println("File: " + file + ", queue: " + context.thread.getClass());
		
		req.status.maybeClearQueue(context.queue);
		context.queue.add(req);
		context.thread.notifyTree();
	}

	public void addResultListener(File file, AreaResultListener listener){
		queues.get(file).listeners.add(listener);
	}

	public void processAreaResult(AreaResult areaResult) {								
		
		for(AreaResultListener listener: queues.get(areaResult.status.file).listeners){
			listener.processAreaResult(areaResult);
		}
	}
}
