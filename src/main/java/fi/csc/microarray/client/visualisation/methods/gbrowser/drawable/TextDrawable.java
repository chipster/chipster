package fi.csc.microarray.client.visualisation.methods.gbrowser.drawable;

import java.awt.Color;

public class TextDrawable extends Drawable {

	private static final int TEXT_HEIGHT = 10;
	public String text;

	public TextDrawable(int x, int y, String text, Color color) {
		super(x, y, color);
		this.text = text;
	}

	public void upsideDown() {
		super.upsideDown();
		y += TEXT_HEIGHT;
	}
}
