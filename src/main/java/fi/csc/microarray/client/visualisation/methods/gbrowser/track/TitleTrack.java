package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;
import java.util.Map;
import java.util.Set;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;

/**
 * Track for placing title texts on top of other tracks.
 *
 */
public class TitleTrack extends Track {

	private Color color;
	private String title;

	public TitleTrack(String title, Color color) {

		this.color = color;
		this.title = title;
		layoutHeight = 10;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();
		drawables.add(new TextDrawable(5, 10, title, color));
		return drawables;
	}

	public void processAreaResult(AreaResult areaResult) {
		// ignore
	}

	@Override
	public int getHeight() {
		return 10;
	}

    @Override
    public Map<AreaRequestHandler, Set<ColumnType>> requestedData() {
        return null;
    }

	@Override
	public String getName() {
		return "title";
	}
}