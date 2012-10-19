package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;
import java.util.Map;
import java.util.Set;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;

/**
 * Empty track can be used for padding.
 *
 */
public class EmptyTrack extends Track {

	public EmptyTrack(View view, int height) {
		super(view, null);
		this.height = height;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();
		drawables.add(new RectDrawable(0, 0, getView().getWidth(), height, Color.WHITE, Color.WHITE));
		return drawables;
	}

	public void processAreaResult(AreaResult areaResult) {
		// ignored
	}
	
    @Override
    public Integer getHeight() {
        return height;
    }

	@Override
	public boolean isStretchable () {
		return false;
	}

    @Override
    public Map<DataSource, Set<ColumnType>> requestedData() {
        return null;
    }

	@Override
	public boolean isConcised() {
		return false;
	}
	
	@Override
	public String getName() {
		return "empty";
	}
}
