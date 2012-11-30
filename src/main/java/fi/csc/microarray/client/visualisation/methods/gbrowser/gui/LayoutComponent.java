package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

public interface LayoutComponent {
	public boolean isFixedHeight();
	public int getHeight();
	public void setHeight(int height);
	public int getMinHeight();
	public boolean isVisible();
	public int getCanvasHeight();
}
