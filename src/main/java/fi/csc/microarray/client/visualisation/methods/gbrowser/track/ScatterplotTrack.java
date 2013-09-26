package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.IndexKey;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.ScatterplotValue;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.SelectionText;
import fi.csc.microarray.util.ScaleUtil;

public class ScatterplotTrack extends Track {

	private static final int TOP_MARGIN = 0;

	private TreeMap<IndexKey, ScatterplotPoint> data = new TreeMap<>();

	private Color defaultColor;
	private int height;
	private float minValue;
	private float maxValue;
	private int minSize = 2; //in pixels
	private boolean xAxisVisible = false;

	private DataType column = null;
	private int floatListIndex;

	private float[] scale;
	
	private ScatterplotTrack(Color color, int height, float minValue, float maxValue) {
		super();
		this.defaultColor = color;
		this.height = height;
		
		int stepCount = ScaleUtil.DEFAULT_STEP_COUNT;
		if (height < 40) {
			stepCount = 2;
		}
				
		scale = ScaleUtil.generateScaleValues(minValue, maxValue, stepCount);
		this.minValue = scale[0];
		this.maxValue = scale[stepCount - 1];
	}


	public ScatterplotTrack(Color color, int height, float minValue, float maxValue, DataType column) {
		this(color, height, minValue, maxValue);
		this.column = column;		
	}

	public ScatterplotTrack(Color color, int height, float minValue, float maxValue,
			int floatArrayIndex, int minBpLength, long maxBpLength) {
		
		this(color, height, minValue, maxValue);
		this.floatListIndex = floatArrayIndex;
	}
	
	public void setMinSize(int minSize) {
		this.minSize = minSize;
	}
	
	public void setXAxisVisible(boolean visible) {
		xAxisVisible = visible;
	}
	
	@Override
	public List<Selectable> getSelectables() {
		List<Selectable> items = new LinkedList<>();

		if (data != null) {
			
			items.add(new PassiveItem(getScaleDrawables()));

			Iterator<IndexKey> iter = data.keySet().iterator();
			while (iter.hasNext()) {
				
				ScatterplotPoint point = data.get(iter.next());

				if (!getView().requestIntersects(point.getRegion())) {
					iter.remove();
					continue;
				}

				point.render(getView(), minSize, this);				 
				items.add(point);
			}
		}

		return items;
	}

	private List<Drawable> getScaleDrawables() {
		
		List<Drawable> drawables = new LinkedList<>();
		
		if (getTrackHeight() > 20) {
			
			//Color scaleColor = new Color(220, 220, 220);
			Color scaleColor = Color.black;
			
			final int CHAR_WIDTH = 6;
			final int START_OFFSET = 3;
			final int END_OFFSET = -13;
			int textOffset = START_OFFSET;
			int maxTextWidth = 0;
			
			for (float scaleValue : scale) {
				
				String text = "" + ScaleUtil.format(scaleValue);
				maxTextWidth = Math.max(maxTextWidth, text.length() * CHAR_WIDTH);					 
			}	
			
			int scaleX = super.getVisibleWidth() - maxTextWidth - 4;
			
			for (float scaleValue : scale) {
				
				int lineY = getY(scaleValue);
				int textY = getY(scaleValue) + 6;
				int lineX2 = scaleX - 6;
				String text = "" + ScaleUtil.format(scaleValue);
				
				if (xAxisVisible && scaleValue == 0f) {
					lineX2 = 0;
				}
				
				drawables.add(new LineDrawable(scaleX, lineY, lineX2, lineY, scaleColor));
				drawables.add(new TextDrawable(scaleX + 4, textY + textOffset, text, scaleColor));
				
				textOffset -= (START_OFFSET - END_OFFSET) / scale.length; 
			}	
			
			drawables.add(new LineDrawable(scaleX, getY(minValue), scaleX, getY(maxValue), scaleColor));		
		}
		return drawables;
	}
	
	public int getY(float value) {
		return (int) ((value - minValue) / (maxValue - minValue) * (height - TOP_MARGIN - 2)  + 1);
	}
	
	public void processDataResult(DataResult dataResult) {
		
		for (Feature feature : dataResult.getFeatures()) {
			
			
			IndexKey key = feature.getIndexKey();
			if (!data.containsKey(key)) {
									
				Float value = null;
				Color itemColor = defaultColor;
				
				SelectionText text = null;
				if (feature.getValueObject() instanceof SelectionText) {
					text = (SelectionText) feature.getValueObject();
				}				
				
				if (column != null) {
					Object obj = feature.values.get(column);
					
					if (obj instanceof SelectionText) {
						SelectionText textObj = (SelectionText) obj;
						text = textObj;
					}
					
					if (obj instanceof ScatterplotValue) {
						ScatterplotValue valueObj = (ScatterplotValue) obj;
						value = valueObj.getScatterplotValue();
						if (valueObj.getScatterplotColor() != null) {
							itemColor = valueObj.getScatterplotColor();
						}

					} else if (obj != null) {
						value = (Float) obj;
					}
				} else {
					
					@SuppressWarnings("unchecked")
					List<Float> floatList = (List<Float>) feature.values.get(DataType.FLOAT_LIST);
					value = floatList.get(floatListIndex);
				}
				
				ScatterplotPoint point = new ScatterplotPoint(feature.region, key, value, itemColor, text);
				if (value != null) {
					data.put(key, point);
				}
			}
		}
	}
	
    @Override
	public void defineDataTypes() {
    	addDataType(DataType.REGION);
    	if (column == null) {
        	addDataType(DataType.FLOAT_LIST);
    	} else {
        	addDataType(column);	
    	}
	}	    
	
	@Override
	public int getTrackHeight() {
		return height;
	}
}
