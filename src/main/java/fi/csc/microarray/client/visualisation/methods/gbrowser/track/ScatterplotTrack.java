package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.IndexKey;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.ScatterplotValue;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.SelectionText;

public class ScatterplotTrack extends ScaleTrack {
	
	private TreeMap<IndexKey, ScatterplotPoint> data = new TreeMap<>();

	private Color defaultColor;
		
	private int minSize = 2; //in pixels

	private DataType column = null;
	private int floatListIndex;

	
	private ScatterplotTrack(Color color, int height, float minValue, float maxValue) {
		super(height, minValue, maxValue);
		this.defaultColor = color;
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
	public void updateLayout() {
		// shared auto scaling makes sense only between read tracks
	}
}
