package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * 
 * Track for showing wig format data
 * 
 */

public class WIGTrack extends Track{
	
	private long maxBpLength;
	private long minBpLength;
	
	private Collection<RegionContent> values = new TreeSet<RegionContent>();
	TreeMap<Long, Float> collector = new TreeMap<Long, Float>();
	private Color color;

	public WIGTrack(View view, DataSource file, Class<? extends AreaRequestHandler> handler, 
			Color color, long minBpLength, long maxBpLength) {
		super(view, file, handler);
		this.color = color;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		
		Collection<Drawable> drawables = getEmptyDrawCollection();
		collector.clear();
		
		Iterator<RegionContent> iter = values.iterator();
		
		// start drawing from the beginning of the track
		float lastX = 1;
		float lastY = 1;
		
        if (values != null) {
        	
        	System.out.println(values.size());
	        while (iter.hasNext()) {
	        	
	            RegionContent value = iter.next();
	
	            // remove those that are not in this view
	            if (!value.region.intercepts(getView().getBpRegion())) {
	                iter.remove();
	                continue;
	            }
	                        	
	            float x = (getView().bpToTrackFloat(value.region.start) +
	                       getView().bpToTrackFloat(value.region.end)) / 2;
	            float y = (Float)value.values.get(ColumnType.VALUE);
	                        
                drawables.add(new LineDrawable((int)lastX, getHeight() - (int)lastY,
                		(int)x, getHeight() - (int)y, color));
                
                lastX = x;
                lastY = y;
	        }
        }

        drawables.add(new LineDrawable((int)lastX, getHeight() - (int)lastY,
        		getView().getWidth(),
        		getHeight()-1, color));
        
		return drawables;
	}
	
	@Override
	public void processAreaResult(AreaResult<RegionContent> areaResult) {
		this.values.add(areaResult.content);
        getView().redraw();
	}

	@Override
	public Collection<ColumnType> getDefaultContents() {
		return Arrays.asList(new ColumnType[] { ColumnType.VALUE });
	}

	@Override
	public boolean isConcised() {
		return false;
	}

	@Override
    public boolean isStretchable() {
        // stretchable unless hidden
        return isVisible();
    }
	
	@Override
	public boolean isVisible(){
		return (super.isVisible() &&
                getView().getBpRegion().getLength() > minBpLength &&
                getView().getBpRegion().getLength() <= maxBpLength);
	}
	
	@Override
    public Integer getHeight() {
        if (isVisible()) {
            return 51;
        } else {
            return 0;
        }
    }

}
