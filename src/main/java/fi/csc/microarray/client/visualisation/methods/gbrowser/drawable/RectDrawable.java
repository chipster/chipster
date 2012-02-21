package fi.csc.microarray.client.visualisation.methods.gbrowser.drawable;

import java.awt.Color;
import java.awt.Rectangle;

import fi.csc.microarray.client.visualisation.methods.gbrowser.TrackContext;

/**
 * A rectangle that can be either filled or not.
 * 
 * @author Petri Klemelä
 *
 */
public class RectDrawable extends Drawable {

	public int width;
	public int height;
	public Color lineColor;
	private String tooltipText;

	public RectDrawable(int x, int y, int width, int height, Color fillColor, Color lineColor, String tooltipText) {
		super(x, y, fillColor);
		this.width = width;
		this.height = height;
		this.lineColor = lineColor;
		this.tooltipText = tooltipText;
	}

	public void upsideDown() {
		super.upsideDown();
		y -= height;
	}

	public RectDrawable(int x, int y, int width, int height, Color fillColor, Color lineColor) {
		this(x, y, width, height, fillColor, lineColor, null);
	}

	public RectDrawable(Rectangle rect, Color fillColor, Color lineColor, String tooltipText) {
		this(rect.x, rect.y, rect.width, rect.height, fillColor, lineColor, tooltipText);
	}

	public RectDrawable(Rectangle rect, Color fillColor, Color lineColor) {
		this(rect.x, rect.y, rect.width, rect.height, fillColor, lineColor, null);
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

	@Override
	public String getTooltipText() {
		return tooltipText;
	}
}
