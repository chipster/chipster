package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.CnaRow;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.CnaRow.Sample;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.IndexKey;

public class CnaFlagTrack extends Track {

	private static final int MIN_VISIBLE_SIZE = 5;

	private static final int SYMBOL_HEIGHT = 10;

	private TreeMap<IndexKey, CnaRow> rows = new TreeMap<IndexKey, CnaRow>();

	private long maxBpLength;
	private long minBpLength;

	private Color gainColor = GBrowserConstants.COLOR_RED;
	private Color lossColor = GBrowserConstants.COLOR_BLUE;

	private int sampleIndex;


	public CnaFlagTrack(Color gainColor, int sampleIndex, Color lossColor, long minBpLength, long maxBpLength) {

		this.sampleIndex = sampleIndex;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		if (rows != null) {
			
			Iterator<IndexKey> iter = rows.keySet().iterator();
			while (iter.hasNext()) {

				CnaRow row = rows.get(iter.next());

				if (!row.getRegion().intersects(getView().getBpRegion())) {
					iter.remove();
					continue;
				}
				
				int y = 0;
				
				Sample sample = row.getSamples().get(sampleIndex);
				
				Color color = Color.lightGray;
				if (sample.getFlag() < 0) {
					color = gainColor;
				} else if (sample.getFlag() > 0) {
					color = lossColor;					
				}

				int alpha = 255 - (int) Math.min(1.0, (Math.abs(sample.getFlag())) * 255.0);

				Color aColor = new Color(color.getRed(), color.getGreen(), color.getBlue(), alpha);

				createDrawable(row.getRegion().start, row.getRegion().end, y, SYMBOL_HEIGHT, aColor, drawables);
			}
		}

		return drawables;
	}

	private void createDrawable(BpCoord startBp, BpCoord endBp, int y, int height, Color c, Collection<Drawable> drawables) {
		Rectangle rect = new Rectangle();

		rect.x = getView().bpToTrack(startBp);
		rect.width = getView().bpToTrack(endBp) - rect.x;		
		
		if (rect.width < MIN_VISIBLE_SIZE) {
			rect.width = MIN_VISIBLE_SIZE;
		}

		rect.y = y;
		rect.height = height;

		drawables.add(new RectDrawable(rect, c, c));
	}

	public void processAreaResult(AreaResult areaResult) {

		for (RegionContent region : areaResult.getContents()) {
			this.rows.put((IndexKey)region.values.get(ColumnType.ID), (CnaRow)region.values.get(ColumnType.VALUE));
		}
		
		getView().redraw();
	}
    
    @Override
    public boolean isVisible() {
        // visible region is not suitable
        return (super.isVisible() &&
                getView().getBpRegion().getLength() > minBpLength &&
                getView().getBpRegion().getLength() <= maxBpLength);
    }
    
    @Override 
    public int getMinHeight() {
    	return SYMBOL_HEIGHT;
    }
}
