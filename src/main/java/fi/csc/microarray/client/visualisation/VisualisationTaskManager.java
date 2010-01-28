package fi.csc.microarray.client.visualisation;

import java.lang.reflect.InvocationTargetException;
import java.util.LinkedList;
import java.util.concurrent.locks.Condition;
import java.util.concurrent.locks.Lock;
import java.util.concurrent.locks.ReentrantLock;

import javax.swing.JComponent;
import javax.swing.SwingUtilities;

import org.apache.log4j.Logger;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.databeans.features.Table;
import fi.csc.microarray.exception.MicroarrayException;
import fi.csc.microarray.util.ThreadUtils;

/**
 * Manages visualisation threads.
 * 
 * @author Petri KlemelÃ¤, Aleksi Kallio
 * 
 */
public class VisualisationTaskManager {

	private class VisualisationRunnable implements Runnable {

		private VisualisationMethodChangedEvent event;

		public VisualisationRunnable(VisualisationMethodChangedEvent event) {
			this.event = event;
		}

		public void run() {
				
				long startTime = System.currentTimeMillis();
				
				// do the actual work (if needed)
				final JComponent visualisation = event.getNewMethod() != null ? frameManager.createVisualisation(event) : null;

				// update GUI in EDT (and wait for it to happen)
				try {
					SwingUtilities.invokeAndWait(new UpdateGuiRunnable(visualisation, event));
				} catch (InterruptedException e) {
					application.reportException(e);
				} catch (InvocationTargetException e) {
					application.reportException(e);
				}
				
/*				logger.debug("Duration, estimated: " + method.estimateDuration(datas) + 
						"\treal: " + (System.currentTimeMillis() - startTime));*/ //NullPointer if method == null
				
				//Enable debug messages to collect stats from visualisation durations
				if(logger.isDebugEnabled() && event.getDatas().size() > 0 && 
						event.getNewMethod() != VisualisationMethod.NONE){ //Possible nullPointers, just for debug
					
					long endTime = System.currentTimeMillis();								
					Table rowCounter = null;
					
					try {
						rowCounter = event.getDatas().get(0).queryFeatures("/column/*").asTable();
					} catch (MicroarrayException e) {
						// Only in debug mode
						e.printStackTrace();
					}
					
					int rowCount = 0;
					
					// 15 000 rows is still reasonably fast, but 500 000 isn't
					while (rowCounter.nextRow()) {
						rowCount++;
					}				

					logger.debug("\tMethod\t" + event.getNewMethod() + 
							"\tDatas\t" + event.getDatas().size() + 
							"\tRows\t" + rowCount + 
							"\tByteLength\t" + event.getDatas().get(0).getContentLength() + 
							"\tTime\t" + (endTime - startTime) + 
							"\tTime/ByteLength\t" + (endTime - startTime)/(float)event.getDatas().get(0).getContentLength() + "\t");
				}
		}
	};

	// update to screen
	private class UpdateGuiRunnable implements Runnable {

		private JComponent visualisation;
		private VisualisationMethodChangedEvent event;

		public UpdateGuiRunnable(JComponent visualisation, VisualisationMethodChangedEvent event) {
			this.visualisation = visualisation;
			this.event = event;
		}

		public void run() {

			frameManager.showVisualisationComponent(visualisation, event);
		}
	};

	private class GuiWorker implements Runnable {

		public void run() {
			if (!Thread.currentThread().isDaemon()) {
				throw new IllegalThreadStateException("GuiWorker must be run on daemon thread");
			}
			
			// this is daemon thread so looping forever does not block JVM exit
			while (true) {
				Runnable visualisationRunnable = null;
				workQueueLock.lock();
				try {
					if (!workQueue.isEmpty()) {
						visualisationRunnable = workQueue.getFirst();
						workQueue.clear(); // pending runnables are cleared
					}
				} finally {
					workQueueLock.unlock();
				}

				if (visualisationRunnable != null) {
					// we had a job to do
					visualisationRunnable.run();
					
				} else {
					// nothing to do, wait for signals
					workQueueLock.lock();
					try {
						guiWorkAvailable.awaitUninterruptibly();
					} finally {
						workQueueLock.unlock();
					}
				}
			}
		}		
	}
	
	private static final Logger logger = Logger.getLogger(VisualisationTaskManager.class);

	private LinkedList<Runnable> workQueue = new LinkedList<Runnable>();
	private Lock workQueueLock = new ReentrantLock();
	private Condition guiWorkAvailable = workQueueLock.newCondition();
	private ClientApplication application = Session.getSession().getApplication();
	private VisualisationFrameManager frameManager;

	public VisualisationTaskManager(VisualisationFrameManager frameManager) {
		this.frameManager = frameManager;
		Thread thread = ThreadUtils.getLowPriorityBackgroundThread(new GuiWorker());
		thread.start();
	}

	public void visualise(VisualisationMethodChangedEvent e) {
		
		workQueueLock.lock();
		try {
			workQueue.add(new VisualisationRunnable(e));
			guiWorkAvailable.signalAll();
			
		} finally {
			workQueueLock.unlock();
		}
	}
}
