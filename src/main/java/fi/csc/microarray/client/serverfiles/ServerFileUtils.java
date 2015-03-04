package fi.csc.microarray.client.serverfiles;

import java.awt.Component;
import java.awt.Container;

import javax.swing.AbstractButton;
import javax.swing.JButton;
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

	
	/**
	 * Find button by tooltip
	 * 
	 * @param parent
	 * @param tooltip
	 */
	private static void hideChildButtonsWithTooltip(Container parent, String tooltip) {
		
		for (Component component : parent.getComponents()) {
			if (component instanceof AbstractButton && tooltip.equals(((AbstractButton)component).getToolTipText())) {
				component.setVisible(false); // hide this
			} else if (component instanceof Container ){
				hideChildButtonsWithTooltip((Container)component, tooltip);
			}
		}
		
	}
	
	/**
	 * Find button by text
	 * 
	 * @param c
	 * @param text
	 * @return
	 */
	private static JButton lookupButton(Container c, String text) {
	    JButton button = null;
	    for (Component comp : c.getComponents()) {
	        if (comp == null) {
	            continue;
	        }
	        if (comp instanceof JButton && (button = (JButton) comp).getText() != null && button.getText().equals(text)) {
	            return button;
	        } else if (comp instanceof Container) {
	            if ((button = lookupButton((Container) comp, text)) != null) {
	                return button;
	            }
	        }
	    }
	    return button;
	}
	
	public static void addButton(JFileChooser fileChooser) {
		// this will initialize the button and its parent component
		fileChooser.setApproveButtonText("approveButton");
		JButton button = lookupButton(fileChooser, "approveButton");
		// reset approve button text
		fileChooser.setApproveButtonText(null);
		// add a new button
		button.getParent().add(new JButton("Save"), 0);
	}
	
	public static void hideApproveButton(JFileChooser fileChooser) {
		// this will initialize the button and its parent component
		fileChooser.setApproveButtonText("approveButton");
		JButton button = lookupButton(fileChooser, "approveButton");
		// reset approve button text
		fileChooser.setApproveButtonText(null);
		// hide it
		button.setVisible(false);
	}
	
	public static void setCancelButtonText(JFileChooser fileChooser, String text) {
		JButton button = lookupButton(fileChooser, "Cancel");
		if (button != null) {
			button.setText(text);
		}
	}
}
