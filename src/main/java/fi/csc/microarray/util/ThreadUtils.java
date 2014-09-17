package fi.csc.microarray.util;

import java.lang.reflect.InvocationTargetException;

import javax.swing.SwingUtilities;

import fi.csc.microarray.client.Session;

/**
 * Factory methods for threads.
 * 
 * @author Aleksi Kallio
 * 
 */
public class ThreadUtils {

	public static Thread getBackgroundThread(Runnable runnable) {
		Thread thread = new Thread(runnable);
		thread.setDaemon(true); // will not block JVM shutdown
		return thread;
	}

	public static Thread getLowPriorityBackgroundThread(Runnable runnable) {
		Thread thread = getBackgroundThread(runnable);
		thread.setPriority(Thread.MIN_PRIORITY + 1); // keep other stuff running
														// smoothly
		return thread;
	}

	/**
	 * Run runnable in EDT and wait for it to complete
	 * 
	 * If this is EDT, call run() directly, otherwise use
	 * SwingUtilities.invokeAndWait().
	 * 
	 * @param runnable
	 */
	public static void runInEDT(Runnable runnable) {
		if (SwingUtilities.isEventDispatchThread()) {
			runnable.run();
		} else {
			try {
				SwingUtilities.invokeAndWait(runnable);
			} catch (InvocationTargetException | InterruptedException e) {
				Session.getSession().getApplication().reportException(e);
			}
		}
	}
}
