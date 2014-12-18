package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import javax.swing.JComponent;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutTool.LayoutMode;

public interface LayoutChild {
	public JComponent getLayoutComponent();
	public LayoutMode getLayoutMode();
	public void updateLayout();
}
