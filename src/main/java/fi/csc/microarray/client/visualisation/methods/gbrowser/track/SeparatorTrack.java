package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;
import java.util.Map;
import java.util.Set;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;

/**
 * Line that separates real tracks.
 *
 */
public class SeparatorTrack extends Track {

	private Color color;
	private long maxBpLength;
	private long minBpLength;
	private int thickness;
	private String name = "separator";
	
	protected SeparatorTrack(GBrowserView view, long minBpLength, long maxBpLength) {
		super(view, null);
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}
	
	public SeparatorTrack(GBrowserView view, Color color, int thickness, long minBpLength, long maxBpLength) {
		this(view, minBpLength, maxBpLength);
		this.color = color;
		this.thickness = thickness;
	}

	public void setName(String name) {
		this.name = name;
	}
	
	@Override
	public String getName() {
		return name;
	}
	
	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();
		for (int i = 0; i < thickness; i++) {
			drawables.add(new LineDrawable(0, i, getView().getWidth(), i, color));
		}

		return drawables;
	}


	public void processAreaResult(AreaResult areaResult) {
		// ignore
	}  

    @Override
    public int getHeight() {
        return thickness;
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
}