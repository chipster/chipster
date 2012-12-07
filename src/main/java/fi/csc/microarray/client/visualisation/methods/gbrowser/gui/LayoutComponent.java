package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutTool.LayoutMode;

public interface LayoutComponent {
	public LayoutMode getLayoutMode();
	public int getHeight();
	public void setHeight(int height);
	public int getMinHeight();
	public boolean isVisible();
	public int getFullHeight();
	public void setLayoutMode(LayoutMode mode);
	public void setDefaultLayoutMode();
}
