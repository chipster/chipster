package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;

/**
 * Empty track can be used for padding.
 *
 */
public class EmptyTrack extends Track {

	public EmptyTrack(int height) {
		super(height);
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();
		drawables.add(new RectDrawable(0, 0, getView().getWidth(), getTrackHeight(), Color.WHITE, Color.WHITE));
		return drawables;
	}

	public void processDataResult(DataResult dataResult) {
		// ignored
	}
	
	@Override
	public String getTrackName() {
		return "empty";
	}
}
