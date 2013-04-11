package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;
import java.util.Map;
import java.util.Set;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;

/**
 * Empty track can be used for padding.
 *
 */
public class EmptyTrack extends Track {

	private Long minBp;

	public EmptyTrack(int height) {

		this.layoutHeight = height;
	}
	
	public EmptyTrack(int height, long minBp) {

		this.layoutHeight = height;
		this.minBp = minBp;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();
		drawables.add(new RectDrawable(0, 0, getView().getWidth(), layoutHeight, Color.WHITE, Color.WHITE));
		return drawables;
	}

	public void processAreaResult(AreaResult areaResult) {
		// ignored
	}

    @Override
    public Map<AreaRequestHandler, Set<ColumnType>> requestedData() {
        return null;
    }
	
	@Override
	public String getName() {
		return "empty";
	}
	
    public boolean isVisible() {
    	if (minBp == null) {
    		return visible;
    	} else {
    		return view.getBpRegion().getLength() >= minBp;
    	}
    }
}
