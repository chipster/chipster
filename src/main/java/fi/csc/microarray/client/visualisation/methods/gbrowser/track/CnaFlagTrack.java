package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.CnaRow;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.CnaRow.Sample;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.IndexKey;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;

public class CnaFlagTrack extends Track {

	private static final int MIN_VISIBLE_SIZE = 5;

	private static final int SYMBOL_HEIGHT = 10;

	private TreeMap<IndexKey, CnaRow> rows = new TreeMap<IndexKey, CnaRow>();

	private Color gainColor = GBrowserConstants.COLOR_RED;
	private Color lossColor = GBrowserConstants.COLOR_BLUE;

	private int sampleIndex;


	public CnaFlagTrack(Color gainColor, int sampleIndex, Color lossColor) {
		super();
		this.sampleIndex = sampleIndex;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		if (rows != null) {
			
			Iterator<IndexKey> iter = rows.keySet().iterator();
			while (iter.hasNext()) {

				CnaRow row = rows.get(iter.next());

				if (!getView().requestIntersects(row.getRegion())) {
					iter.remove();
					continue;
				}
				
				int y = 0;
				
				Sample sample = row.getSamples().get(sampleIndex);
				
				Color color = Color.lightGray;
				if (sample.getFlag() != null) {
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

	public void processDataResult(DataResult dataResult) {

		for (Feature region : dataResult.getFeatures()) {
			this.rows.put((IndexKey)region.values.get(DataType.ID), (CnaRow)region.values.get(DataType.VALUE));
		}
	}
	
	@Override
	public void defineDataTypes() {
		addDataType(DataType.ID);
		addDataType(DataType.VALUE);
	}    
    
    @Override 
    public int getTrackHeight() {
    	return SYMBOL_HEIGHT;
    }
}
