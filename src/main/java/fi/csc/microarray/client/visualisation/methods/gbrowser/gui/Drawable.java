package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.Color;
import java.awt.Graphics2D;

import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackContext;

/**
 * Base class for all drawable objects (shapes). Actual drawing is done by {@link GBrowserView} objects. 
 * 
 * @author Petri Klemelä
 *
 */
public abstract class Drawable {

	public int x;
	public int y;
	public Color color;

	public Drawable(int x, int y, Color color) {
		this.x = x;
		this.y = y;
		this.color = color;
	}
	
	public Drawable(Color color) {
		this.x = -1;
		this.y = -1;
		this.color = color;
	}
	
	public void draw(Graphics2D g, int x, int y) {

		g.setPaint(this.color);
	}	

	public void upsideDown() {
		y = -y;
	}
	
	/**
     * @return minimum y value of all pixels drawn.
     */
    public int getMaxY() {
        return y;
    }
    
    /**
     * Expand this drawable to fit into given context.
     */
    public void expand(TrackContext context) { }
}
