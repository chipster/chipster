package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Generic track for showing high level distribution of items (genes, transcripts, reads...) on the genome.
 * The result is an approximation.
 *
 */
public class TabixIntensityTrack extends Track {

	private SortedSet<RegionContent> values = new TreeSet<RegionContent>();
	private long minBpLength;
	private long maxBpLength;
	private Color color;

	public TabixIntensityTrack(View view, DataSource file, Class<? extends AreaRequestHandler> handler,
	        Color c, long minBpLength, long maxBpLength) {
		super(view, file, handler);
		this.color = c;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {

		Collection<Drawable> drawables = getEmptyDrawCollection();

		Iterator<RegionContent> iterator = values.iterator();
		
		int lastX = 0;
		int lastY = (int) getView().getTrackHeight();
		
		while (iterator.hasNext()) {

			RegionContent regCont = iterator.next();
			
			// remove values that have gone out of view
			if (!regCont.region.intercepts(getView().getBpRegion())) {
				iterator.remove();
				continue;
			}
			
			// do the plotting for this consised value
			int x1 = getView().bpToTrack(regCont.region.start);
			//int x2 = getView().bpToTrack(regCont.region.end) + 2;
			int y2 = (int) getView().getTrackHeight();						
			
			int val = (int) Math.min( (Math.log((Float)regCont.values.get(ColumnType.VALUE)) / regCont.region.getLength() * 1000), 
					getView().getTrackHeight());
			
			int y1 = (int) (-val + y2);
			
			drawables.add(new LineDrawable(lastX, lastY, x1, y1, color));
			
			lastX = x1;
			lastY = y1;
		}

		return drawables;
	}

	public void processAreaResult(AreaResult<RegionContent> areaResult) {
		
		
		if (areaResult.status.concise == this.isConcised() && 
				areaResult.content.values.get(ColumnType.VALUE) != null &&
				areaResult.content.values.get(ColumnType.VALUE) instanceof Float &&
				areaResult.content.region.intercepts(getView().getBpRegion())) { 
			
			values.add(areaResult.content);
			getView().redraw();
		}
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
