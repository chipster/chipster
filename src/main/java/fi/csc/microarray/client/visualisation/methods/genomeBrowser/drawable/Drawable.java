package fi.csc.microarray.client.visualisation.methods.genomeBrowser.drawable;

import java.awt.Color;

public class Drawable {

	public int x;
	public int y;
	public Color color;
	
	public Drawable(int x, int y, Color color){
		this.x = x;
		this.y = y;
		this.color = color;
	}

	public void upsideDown() {
		y = -y;
	}
}
