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
    public int getMinY() {
        return y;
    }
    
    @Override
    public void expand(TrackContext context) {
        // firstly coordinates are converted from user coordinate space to
        // euclidian space, then multiplied by expansion ratio and converted back
        int maxY = context.trackHeight-1;
        this.y = Math.max(maxY - Math.round((-this.y + maxY) *
                context.expansionRatio), 0);
        
//        int y2 = this.y + this.height;
//        
//        int newY2 = Math.max(maxY - Math.round((-y2 + maxY) *
//                context.expansionRatio), 0);
//        
//        this.height =  newY2 - y;
        
        int newHeight = height = (int) (this.height * context.expansionRatio);
        y -= newHeight-height;
        this.height = newHeight;
    }
}
