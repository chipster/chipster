package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.stack.IndexKey;

public class ScatterplotTrack extends Track {

	private static final int TOP_MARGIN = 0;

	private TreeMap<IndexKey, RegionContent> data = new TreeMap<IndexKey, RegionContent>();

	private long maxBpLength;
	private long minBpLength;

	private Color color;
	private int height;
	private float minValue;
	private float maxValue;

	private ColumnType column = null;
	private int floatListIndex;


	public ScatterplotTrack(Color color, int height, float minValue, float maxValue, ColumnType column, long minBpLength, long maxBpLength) {

		this.color = color;
		this.height = height;
		this.minValue = minValue;
		this.maxValue = maxValue;
		this.column = column;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;		
	}

	public ScatterplotTrack(Color color, int height, float minValue, float maxValue,
			int floatArrayIndex, int minBpLength, long maxBpLength) {
		
		this.color = color;
		this.height = height;
		this.minValue = minValue;
		this.maxValue = maxValue;
		this.floatListIndex = floatArrayIndex;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		if (data != null) {
			
			if (getMinHeight() > 20) {
				drawables.add(getRectDrawable(0, 17, minValue, Color.gray));
				drawables.add(getRectDrawable(0, 17, maxValue, Color.gray));
				drawables.add(new TextDrawable(1, 12, "" + minValue, Color.gray));
				drawables.add(new TextDrawable(1, height - TOP_MARGIN - 2, "" + maxValue, Color.gray));
			}

			Iterator<IndexKey> iter = data.keySet().iterator();
			while (iter.hasNext()) {

				RegionContent regionContent = data.get(iter.next());

				if (!regionContent.region.intersects(getView().getBpRegion())) {
					iter.remove();
					continue;
				}
				
				
				int x1 = getView().bpToTrack(regionContent.region.start);
				int width = getView().bpToTrack(regionContent.region.end) - x1;
				
				width = Math.max(width, 2);
				
				float value;
				
				if (column != null) {
					value = (Float) regionContent.values.get(column);
				} else {
					
					@SuppressWarnings("unchecked")
					List<Float> floatList = (List<Float>) regionContent.values.get(ColumnType.FLOAT_LIST);
					value = floatList.get(floatListIndex);
				}
				
				if (value >= minValue && value <= maxValue) {
					drawables.add(getRectDrawable(x1, width, value, color));
				}
			}
		}

		return drawables;
	}
	
	public RectDrawable getRectDrawable(int x1, int width, float value, Color color) {
		
		int y = (int) ((value - minValue) / (maxValue - minValue) * (height - TOP_MARGIN - 2));
								
		Rectangle rect = new Rectangle();

		rect.y = y;
		rect.x = x1;
		rect.width = width;
		rect.height = 2;

		return new RectDrawable(rect, color, color);
	}
	
	public void processAreaResult(AreaResult areaResult) {

		for (RegionContent region : areaResult.getContents()) {
			this.data.put((IndexKey)region.values.get(ColumnType.ID), region);
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
		return height;
	}
}
