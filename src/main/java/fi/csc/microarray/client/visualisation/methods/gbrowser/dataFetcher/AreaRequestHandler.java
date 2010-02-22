package fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher;
import java.util.Queue;

import javax.swing.SwingUtilities;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequest;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;


public abstract class AreaRequestHandler extends Thread{
		
		private Queue<AreaRequest> areaRequestQueue;
		private AreaResultListener areaResultListener;
		
		public AreaRequestHandler(Queue<AreaRequest> areaRequestQueue, 
				AreaResultListener areaResultListener) {
			
			super();		
			this.areaRequestQueue = areaRequestQueue;		
			this.areaResultListener = areaResultListener;
			this.setDaemon(true);
		}

		public synchronized void run(){
			
			while(true){
				AreaRequest areaRequest;
				if((areaRequest = areaRequestQueue.poll()) != null){
					areaRequest.status.areaRequestCount = areaRequestQueue.size();
					processAreaRequest(areaRequest);
				}
				
				boolean isWorkToDo = checkOtherQueues();
				
				if(areaRequest == null && !isWorkToDo){
					try {
						this.wait();
					} catch (InterruptedException e) {
						//Just poll the queues and wait again
					}
				}
			}
		}
		
		/**
		 * @return true if should be called again before putting thread to wait
		 */
		protected boolean checkOtherQueues() {
			//Hook for checking fileFetcherQueue
			return false;
		}

		public synchronized void notifyTree(){
			notifyAll();
		}

		protected abstract void processAreaRequest(AreaRequest areaRequest);
		
		public void createAreaResult(final AreaResult areaResult) {
			
			SwingUtilities.invokeLater(new Runnable(){
				public void run(){						
					areaResultListener.processAreaResult(areaResult);
				}
			});
		}
}
