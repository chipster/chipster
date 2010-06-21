package fi.csc.microarray.client.visualisation.methods.gbrowser.drawable;

import java.awt.Color;

import fi.csc.microarray.client.visualisation.methods.gbrowser.TrackContext;

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

	public String toString() {
		return "LineDrawable (" + x + ", " + y + ") - (" + x2 + ", " + y2 + ")";
	}

    @Override
    public int getMinY() {
        return Math.min(y, y2);
    }
    
    @Override
    public void expand(TrackContext context) {
        // firstly coordinates are converted from user coordinate space to
        // euclidian space, then multiplied by expansion ratio and converted back
        this.y = context.trackHeight - Math.round((-this.y + context.trackHeight) *
                context.expansionRatio);
        this.y2 = context.trackHeight - Math.round((-this.y2 + context.trackHeight) *
                context.expansionRatio);
    }
}
