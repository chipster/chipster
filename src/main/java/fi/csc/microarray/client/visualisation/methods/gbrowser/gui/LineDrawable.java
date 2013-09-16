package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.Color;
import java.awt.Graphics2D;

import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackContext;

/**
 * A line that goes from one point to another.
 * 
 * @author Petri Klemelä
 *
 */
public class LineDrawable extends Drawable {

	public int x2;
	public int y2;

	public LineDrawable(int x, int y, int x2, int y2, Color color) {
		super(x, y, color);
		this.x2 = x2;
		this.y2 = y2;
	}
	
	public void draw(Graphics2D g, int x, int y) {
		super.draw(g, x, y);
		g.drawLine(this.x + x, this.y + y, this.x2 + x, this.y2 + y);
	}

	public void upsideDown() {
		super.upsideDown();
		y2 = -y2;
	}

	public String toString() {
		return "LineDrawable (" + x + ", " + y + ") - (" + x2 + ", " + y2 + ")";
	}

    @Override
    public int getMaxY() {
        return Math.min(y, y2);
    }
    
    @Override
    public void expand(TrackContext context) {

        int maxY = context.trackHeight-1;
        this.y = Math.min(maxY, (int) (this.y * context.expansionRatio));
        this.y2 = Math.min(maxY, (int) (this.y2 * context.expansionRatio));
    }
}
