package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;

/**
 * Line that separates real tracks.
 *
 */
public class SeparatorTrack extends Track {

	private Color color;
	private int thickness;
	private String name = "separator";
	
	public SeparatorTrack(Color color, int thickness) {
		super();
		this.color = color;
		this.thickness = thickness;
	}

	public void setTrackName(String name) {
		this.name = name;
	}
	
	@Override
	public String getTrackName() {
		return name;
	}
	
	@Override
	public Collection<Drawable> getDrawables() {
		
		Collection<Drawable> drawables = getEmptyDrawCollection();
		for (int i = 0; i < thickness; i++) {
			drawables.add(new LineDrawable(0, i + 1, getView().getWidth(), i + 1, color));		
		}

		return drawables;
	}

	public void processDataResult(DataResult dataResult) {
		// ignore
	}  

    @Override
    public int getTrackHeight() {
        return thickness;
    }    
}