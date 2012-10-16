package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
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
	private Color color;

	public WIGTrack(View view, DataSource file,	Color color, long minBpLength, long maxBpLength) {
		super(view, file);
		this.color = color;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		
		Collection<Drawable> drawables = getEmptyDrawCollection();
		
		Iterator<RegionContent> iter = values.iterator();
		
		// start drawing from the beginning of the track
		float lastX = 1;
		float lastY = 1;
		
        if (values != null) {
        	
	        while (iter.hasNext()) {
	        	
	            RegionContent value = iter.next();
	
	            // remove those that are not in this view
	            if (!value.region.intersects(getView().getBpRegion())) {
	                iter.remove();
	                continue;
	            }
	                        	
	            float x = (getView().bpToTrackFloat(value.region.start) +
	                       getView().bpToTrackFloat(value.region.end)) / 2;
	            float y = Float.parseFloat((String)value.values.get(ColumnType.VALUE));
	                        
                drawables.add(new LineDrawable((int)lastX, getHeight() - (int)lastY,
                		(int)x, getHeight() - (int)y, color));
                
                lastX = x;
                lastY = y;
	        }
        }

        // finish the line till the end of the screen
        drawables.add(new LineDrawable((int)lastX, getHeight() - (int)lastY,
        		getView().getWidth(),
        		getHeight()-1, color));
        
		return drawables;
	}
	
	@Override
	public void processAreaResult(AreaResult areaResult) {
		this.values.addAll(areaResult.getContents());
        getView().redraw();
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

	@Override
	public Map<DataSource, Set<ColumnType>> requestedData() {
		HashMap<DataSource, Set<ColumnType>> datas = new
        HashMap<DataSource, Set<ColumnType>>();
        datas.put(file, new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {
                ColumnType.VALUE })));
        return datas;
	}

}
