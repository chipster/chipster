package fi.csc.microarray.util;

import javax.swing.SwingUtilities;

public class SwingTools {

	public static void checkEDT() {
		if (!SwingUtilities.isEventDispatchThread()) {
			throw new RuntimeException("outside EDT");
		}
	}
	
	
	public static void runInEventDispatchThread(Runnable runnable) {
		try {
			if (SwingUtilities.isEventDispatchThread()) {
				runnable.run();
				
			} else {
				SwingUtilities.invokeAndWait(runnable);
			}

		} catch (Exception e) {
			throw new RuntimeException(e);
		}
	}
}
