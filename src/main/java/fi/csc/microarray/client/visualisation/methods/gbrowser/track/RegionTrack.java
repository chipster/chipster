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
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.SelectionText;

/**
 * Track for showing the location of predicted peaks. Peaks cannot overlap. 
 *
 */
public class RegionTrack extends Track {

	private TreeMap<IndexKey, PeakSelectable> peaks = new TreeMap<>();
	
	private Color color;

	public RegionTrack(Color color) {
		super();
		this.color = color;
	}
	
	@Override
	public List<Selectable> getSelectables() {
		List<Selectable> items = new LinkedList<>();

		if (peaks != null) {

			Iterator<IndexKey> iter = peaks.keySet().iterator();
			while (iter.hasNext()) {
				
				PeakSelectable selectable = peaks.get(iter.next());

				if (!getView().requestIntersects(selectable.getRegion())) {
					iter.remove();
					continue;
				}

				selectable.render(getView(), color);				 
				items.add(selectable);
			}
		}

		return items;
	}

	public void processDataResult(DataResult dataResult) {
		
		for (Feature feature : dataResult.getFeatures()) {
			
			IndexKey key = feature.getIndexKey();

			if (!peaks.containsKey(key)) {
				SelectionText text = null;

				Object value = feature.getValueObject();
				if (value instanceof SelectionText) {
					text = (SelectionText) value;
				}
				PeakSelectable selectable = new PeakSelectable(feature.region, key, text);
				peaks.put(key, selectable);
			}
		}
	}    
	
	@Override
	public void defineDataTypes() {
		addDataType(DataType.CHROMOSOME);
		addDataType(DataType.START);
		addDataType(DataType.END);
	}
	
	@Override
	public int getTrackHeight() {
		return PeakSelectable.PEAK_SYMBOL_HEIGHT + 1;
	}
}
