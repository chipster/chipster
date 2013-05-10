package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;

/**
 * Track for placing title texts on top of other tracks.
 *
 */
public class TitleTrack extends Track {

	private Color color;
	private String title;
	private Color bgColor;

	public TitleTrack(String title, Color color) {

		this.color = color;
		this.title = title;
		layoutHeight = 10;
	}
	
	public TitleTrack(String title, Color color, Color bgColor) {

		this(title, color);
		this.bgColor = bgColor;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();
		
		if (bgColor != null) {
			drawables.add(new RectDrawable(0, 0, view.getWidth(), getHeight(), bgColor, bgColor));
		}
		
		drawables.add(new TextDrawable(5, 10, title, color));				
		return drawables;
	}

	public void processDataResult(DataResult dataResult) {
		// ignore
	}

	@Override
	public int getHeight() {
		return 10;
	}

	@Override
	public String getName() {
		return "title";
	}
}