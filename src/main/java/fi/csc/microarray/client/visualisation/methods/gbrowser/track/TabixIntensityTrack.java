package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeMap;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Generic track for showing high level distribution of items (genes, transcripts, reads...) on the genome.
 */
public class TabixIntensityTrack extends Track {

	private SortedSet<RegionContent> values = new TreeSet<RegionContent>();
	private long minBpLength;
	private long maxBpLength;
	private Color color;

	public TabixIntensityTrack(View view, DataSource file, Color c, long minBpLength, long maxBpLength) {
		super(view, file);
		this.color = c;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}
	
	@Override
	public Collection<Drawable> getDrawables() {

		Collection<Drawable> drawables = getEmptyDrawCollection();

		Iterator<RegionContent> iterator = values.iterator();
		
		
		// Construct coverage profile by finding where it changes (derivate)
		TreeMap<Integer, Integer> profileDerivate = new TreeMap<Integer, Integer>();
		while (iterator.hasNext()) {

			RegionContent regCont = iterator.next();
			
			// Remove values that have gone out of view
			if (!regCont.region.intersects(getView().getBpRegion())) {
				iterator.remove();
				continue;
			}
			
			// Process this read
			int startX = getView().bpToTrack(regCont.region.start);
			int endX = getView().bpToTrack(regCont.region.end);
			int height = (int) Math.min( (Math.log((Float)regCont.values.get(ColumnType.VALUE)) / regCont.region.getLength() * 1000), 
					getView().getTrackHeight());
			
			// Update profile at the start location
			if (profileDerivate.get(startX) == null) {
				profileDerivate.put(startX, 0);
			}
			profileDerivate.put(startX, profileDerivate.get(startX) + height);
			
			// Update profile at the end location
			if (profileDerivate.get(endX) == null) {
				profileDerivate.put(endX, 0);
			}
			profileDerivate.put(endX, profileDerivate.get(endX) - height);
		}
		

		// Iterate over the profile changes (derivate) and draw the actual profile
		Iterator<Integer> profileIterator = profileDerivate.keySet().iterator();
		int lastX = -1;
		int lastY = 0;
		int lastDrawnX = -1;
		int lastDrawnY = -1;
		while (profileIterator.hasNext()) {

			Integer x = profileIterator.next();
			Integer newY = lastY + profileDerivate.get(x);

			// Smooth profile by removing less than 10 pixel deviations
			if (lastX != -1 && (x-lastX) > 10) {
				
				if (lastDrawnX != -1) {
					drawables.add(new LineDrawable(lastDrawnX, lastDrawnY, lastX, newY, color));
				}
				
				drawables.add(new LineDrawable(lastX, newY, x, newY, color));
				lastDrawnX = x;
				lastDrawnY = newY;
			}

			lastX = x;
			lastY = newY;
		}

		return drawables;
	}

	public void processAreaResult(AreaResult areaResult) {
		
		for (RegionContent content : areaResult.getContents()) {
			if (areaResult.getStatus().concise == this.isConcised() && 
					content.values.get(ColumnType.VALUE) != null &&
					content.values.get(ColumnType.VALUE) instanceof Float &&
					content.region.intersects(getView().getBpRegion())) { 

				values.add(content);
			}
		}
		
		getView().redraw();
	}

	public Integer getHeight() {
		if (isVisible()) {
			return super.getHeight();
		} else {
			return 0;
		}
	}
	   
    @Override
    public boolean isStretchable() {
        // stretchable unless hidden
        return isVisible();
    }
    
    @Override
    public boolean isVisible() {
        // visible region is not suitable
        return (super.isVisible() &&
                getView().getBpRegion().getLength() > minBpLength &&
                getView().getBpRegion().getLength() < maxBpLength);
    }

    @Override
    public Map<DataSource, Set<ColumnType>> requestedData() {
        HashMap<DataSource, Set<ColumnType>> datas = new
                HashMap<DataSource, Set<ColumnType>>();
        datas.put(file, new HashSet<ColumnType>());
        return datas;
    }

	@Override
	public boolean isConcised() {
		return false;
	}
}
