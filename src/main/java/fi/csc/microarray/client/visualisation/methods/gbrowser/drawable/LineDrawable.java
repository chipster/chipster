package fi.csc.microarray.client.visualisation.methods.gbrowser.drawable;

import java.awt.Color;


public class LineDrawable extends Drawable {

	public int x2;
	public int y2;
	
	public LineDrawable(int x, int y, int x2, int y2, Color color) {
		super(x, y, color);
		this.x2 = x2;
		this.y2 = y2;
	}
	
	public void upsideDown() {
		super.upsideDown();
		y2 = -y2;
	}
	
	public String toString(){
		return "LineDrawable (" + x + ", " +  y + ") - (" + x2 + ", " + y2 + ")";
	}
}
