package fi.csc.microarray.client.serverfiles;

import java.awt.Component;
import java.awt.Container;

import javax.swing.AbstractButton;
import javax.swing.JFileChooser;

public class ServerFileUtils {

	
	public static void hideJFileChooserButtons(JFileChooser sessionFileChooser) {
		// stupid graphical buttons do not seem to have anything better than tooltip for identification
		hideChildButtonsWithTooltip(sessionFileChooser, "Sessions at server");
		hideChildButtonsWithTooltip(sessionFileChooser, "Up One Level");
		hideChildButtonsWithTooltip(sessionFileChooser, "Remote sessions");
		hideChildButtonsWithTooltip(sessionFileChooser, "Create New Folder");
		hideChildButtonsWithTooltip(sessionFileChooser, "List");
		hideChildButtonsWithTooltip(sessionFileChooser, "Details");
	}

	
	private static void hideChildButtonsWithTooltip(Container parent, String tooltip) {
		
		for (Component component : parent.getComponents()) {
			if (component instanceof AbstractButton && tooltip.equals(((AbstractButton)component).getToolTipText())) {
				component.setVisible(false); // hide this
			} else if (component instanceof Container ){
				hideChildButtonsWithTooltip((Container)component, tooltip);
			}
		}
		
	}

}
