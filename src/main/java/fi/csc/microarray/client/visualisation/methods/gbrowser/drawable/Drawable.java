package fi.csc.microarray.client.visualisation.methods.gbrowser.drawable;

import java.awt.Color;

public class Drawable {

	public int x;
	public int y;
	public Color color;

	public Drawable(int x, int y, Color color) {
		this.x = x;
		this.y = y;
		this.color = color;
	}

	public void upsideDown() {
		y = -y;
	}
	
	/**
	 * @return maximum y value of all pixels drawn.
	 */
	public int getMaxY() {
	    return y;
	}
}
