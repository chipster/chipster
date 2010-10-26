package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;
import java.util.Map;
import java.util.Set;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Line that separates real tracks.
 *
 */
public class SeparatorTrack extends Track {

	private Color color;
	private boolean thick;
	private long maxBpLength;
	private long minBpLength;

	public SeparatorTrack(View view, Color color, boolean thick, long minBpLength, long maxBpLength) {
		super(view, null);
		this.color = color;
		this.thick = thick;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();
		if (thick) {
			drawables.add(new LineDrawable(0, 1, getView().getWidth(), 1, new Color(243, 243, 243)));
			drawables.add(new LineDrawable(0, 2, getView().getWidth(), 2, new Color(224, 224, 224)));
			drawables.add(new LineDrawable(0, 3, getView().getWidth(), 3, new Color(176, 176, 176)));
			drawables.add(new LineDrawable(0, 4, getView().getWidth(), 4, new Color(64, 64, 64)));
			// White space here
			drawables.add(new LineDrawable(0, 7, getView().getWidth(), 7, new Color(64, 64, 64)));
			drawables.add(new LineDrawable(0, 8, getView().getWidth(), 8, new Color(176, 176, 176)));
			drawables.add(new LineDrawable(0, 9, getView().getWidth(), 9, new Color(224, 224, 224)));
			drawables.add(new LineDrawable(0, 10, getView().getWidth(), 10, new Color(243, 243, 243)));
			
		} else {
			drawables.add(new LineDrawable(0, 0, getView().getWidth(), 0, color));
		}

		return drawables;
	}

	public void processAreaResult(AreaResult<RegionContent> areaResult) {
		// ignore
	}  

    @Override
    public Integer getHeight() {
        if (isVisible()) {
            return thick ? 11 : 1;
            
        } else {
            return 0;
        }
    }
    
    @Override
    public boolean isStretchable() {
        return false;
    }
    
    @Override
    public boolean isVisible() {
        // visible region is not suitable
        return (super.isVisible() &&
                getView().getBpRegion().getLength() > minBpLength &&
                getView().getBpRegion().getLength() <= maxBpLength);
    }    

    @Override
    public Map<DataSource, Set<ColumnType>> requestedData() {
        return null;
    }

	@Override
	public boolean isConcised() {
		return false;
	}
}