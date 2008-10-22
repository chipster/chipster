package fi.csc.microarray.util;

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
		thread.setPriority(Thread.MIN_PRIORITY + 1); // keep other stuff running smoothly
		return thread;
	}

}
