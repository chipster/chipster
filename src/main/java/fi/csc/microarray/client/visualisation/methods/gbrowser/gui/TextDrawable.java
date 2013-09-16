package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.awt.Color;
import java.awt.Graphics2D;

/**
 * A single row of text.
 * 
 * @author Petri Klemelä
 *
 */
public class TextDrawable extends Drawable {

	private static final int TEXT_HEIGHT = 8;
	public String text;

	public TextDrawable(int x, int y, String text, Color color) {
		super(x, y, color);
		this.text = text;
	}
	
	public void draw(Graphics2D g, int x, int y) {

		super.draw(g, x, y);

		g.setFont(g.getFont().deriveFont(10f));		
		g.drawString(this.text, this.x + x, this.y + y);
	}

	public void upsideDown() {
		super.upsideDown();
		y += TEXT_HEIGHT;
	}
    
    @Override
    public int getMaxY() {
        return y - TEXT_HEIGHT;
    }
}
