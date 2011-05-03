package fi.csc.microarray.util;

import java.awt.Color;

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
	
	public static String colorToHexString(Color color) {
        String colorString = Integer.toHexString(color.getRGB());
        return "#" + colorString.substring(2, colorString.length());
	}

	public static Color hexStringToColor(String hexString) {
        return Color.decode(hexString);
	}


}
