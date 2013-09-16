package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;

import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackContext;

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

	public RectDrawable(int x, int y, int width, int height, Color fillColor, Color lineColor, String tooltipText) {
		super(x, y, fillColor);
		this.width = width;
		this.height = height;
		this.lineColor = lineColor;
	}
	
	public RectDrawable(int x, int y, int width, int height, Color fillColor, Color lineColor) {
		this(x, y, width, height, fillColor, lineColor, null);
	}

	public RectDrawable(Rectangle rect, Color fillColor, Color lineColor) {
		this(rect.x, rect.y, rect.width, rect.height, fillColor, lineColor, null);
	}
	
	public void draw(Graphics2D g, int x, int y) {
		
		super.draw(g, x, y);

		// Draw fill
		if (this.color != null) {
			g.setPaint(this.color);
			g.fillRect(this.x + x, this.y + y, this.width, this.height);
		}

		// Draw outline after fill to hide gaps between adjacent rectangles
		if (this.lineColor != null) {
			g.setColor(this.lineColor);
			g.drawRect(this.x + x, this.y + y, this.width - 1, this.height - 1);
		}
	}


	public void upsideDown() {
		super.upsideDown();
		y -= height;
	}
    
    @Override
    public int getMaxY() {
        return y + height;
    }
    
    @Override
    public void expand(TrackContext context) {

        int maxY = context.trackHeight - 1;
        int newY = (int) (this.y * context.expansionRatio);
        int newHeight = (int) (this.height * context.expansionRatio);
        
        this.y = Math.min(newY, maxY);
        this.height = Math.min(newHeight, maxY - y);
    }

	public Rectangle getBounds() {
		return new Rectangle(x, y, width, height);
	}
}
