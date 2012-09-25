package fi.csc.microarray.client.visualisation;

import java.lang.reflect.InvocationTargetException;
import java.util.Iterator;
import java.util.LinkedList;

import javax.swing.JComponent;
import javax.swing.SwingUtilities;

import fi.csc.microarray.client.ClientApplication;
import fi.csc.microarray.client.Session;
import fi.csc.microarray.util.ThreadUtils;

/**
 * Manages visualisation threads.
 * 
 * @author Petri Klemelä, Aleksi Kallio
 * 
 */
public class VisualisationTaskManager {

	private class VisualisationRunnable implements Runnable {

		private VisualisationMethodChangedEvent event;
		private boolean abandonedThread = false;

		public VisualisationRunnable(VisualisationMethodChangedEvent event) {
			this.event = event;
		}
		
		public void abandonThread() {
			abandonedThread = true;
		}

		public void run() {

			// do the actual work (if needed)
			try {
				final JComponent visualisation = event.getNewMethod() != null ? frameManager.createVisualisation(event) : null;

				// update GUI in EDT (and wait for it to happen)
				if (!abandonedThread) {
					SwingUtilities.invokeAndWait(new UpdateGuiRunnable(visualisation, event));
				}
			} catch (InterruptedException e) {
				application.reportException(e);
			} catch (InvocationTargetException e) {
				application.reportException(e);
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
	

	private LinkedList<VisualisationRunnable> visualisationRunnables = new LinkedList<VisualisationRunnable>();

	private ClientApplication application = Session.getSession().getApplication();
	private VisualisationFrameManager frameManager;

	public VisualisationTaskManager(VisualisationFrameManager frameManager) {
		this.frameManager = frameManager;
	}

	public void visualise(VisualisationMethodChangedEvent e) {
		
		/**
		 * When a new visualisation is requested, the old visualisation thread may still be running. Because it isn't
		 * possible to stop those old visualisation threads, we just flag them abandoned. If the old thread finishes later,
		 * this flag prevents it form updating the gui, where the user is already using the new visualisation. If
		 * the thread is really stuck, there isn't much we can do, it will be ended when the application is closed.
		 * 
		 * This list visualisationRunnables contains the threads (or actually the runnables that run in the threads) that 
		 * have been started, but aren't yet flagged abandoned. Usually the thread should have finished already when 
		 * the next visualisation is started, but those threads are flagged anyway. Flagged threads are removed from the list.
		 */
		Iterator<VisualisationRunnable> iter = visualisationRunnables.iterator();	
		while (iter.hasNext()) {
			
			VisualisationRunnable vr = iter.next();
			vr.abandonThread();
			iter.remove();
		}
		
		VisualisationRunnable runnable = new VisualisationRunnable(e);
		visualisationRunnables.add(runnable);
		
		Thread thread = ThreadUtils.getLowPriorityBackgroundThread(runnable);
		thread.start();
	}
}
