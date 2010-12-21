package fi.csc.microarray.client.visualisation.methods.gbrowser.drawable;

import java.awt.Color;
import java.awt.Rectangle;

import fi.csc.microarray.client.visualisation.methods.gbrowser.TrackContext;

public class RectDrawable extends Drawable {

	public int width;
	public int height;
	public Color lineColor;

	public RectDrawable(int x, int y, int width, int height, Color fillColor, Color lineColor) {
		super(x, y, fillColor);
		this.width = width;
		this.height = height;
		this.lineColor = lineColor;
	}

	public void upsideDown() {
		super.upsideDown();
		y -= height;
	}

	public RectDrawable(Rectangle rect, Color fillColor, Color lineColor) {
		this(rect.x, rect.y, rect.width, rect.height, fillColor, lineColor);
	}
    
    @Override
    public int getMaxY() {
        return y + height;
    }
    
    @Override
    public void expand(TrackContext context) {

        int maxY = context.trackHeight-1;
        
        int newHeight = height = (int) (this.height * context.expansionRatio);
        
        this.height = Math.min(newHeight, maxY);
    }
}
