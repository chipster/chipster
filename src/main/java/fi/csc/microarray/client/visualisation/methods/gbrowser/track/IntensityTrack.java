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
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Generic track for showing high level distribution of items (genes, transcripts, reads...) on the genome.
 * The result is an approximation.
 *
 */
public class IntensityTrack extends Track {

	private SortedSet<RegionContent> values = new TreeSet<RegionContent>();
	private long minBpLength;
	private Color color;

	public IntensityTrack(View view, DataSource file, Class<? extends AreaRequestHandler> handler,
	        Color c, long maxBpLength) {
		super(view, file, handler);
		this.color = c;
		this.minBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {

		Collection<Drawable> drawables = getEmptyDrawCollection();

		Iterator<RegionContent> iterator = values.iterator();
		while (iterator.hasNext()) {

			RegionContent regCont = iterator.next();
			
			// remove values that have gone out of view
			if (!regCont.region.intercepts(getView().getBpRegion())) {
				iterator.remove();
				continue;
			}
			
			// do the plotting for this consised value
			int x1 = getView().bpToTrack(regCont.region.start);
			int x2 = getView().bpToTrack(regCont.region.end) + 2;
			int y2 = (int) getView().getTrackHeight();						
			
			int val = (int) Math.min(Math.log((Float) (regCont.values.get(ColumnType.VALUE))) * 4, getView().getTrackHeight() / 4);
			int y1 = (int) (-val + y2);

			// FIXME implement overlaying tracks; currently is drawn above the track,
			//       so it would merge with the previous track
			// FIXME when region content value is close to zero, y1 can be something
			//       like -2147483605, which we probably don't want to plot
			drawables.add(new RectDrawable(x1, y1, x2 - x1, y2 - y1, color, null));

		}

		return drawables;
	}

	public void processAreaResult(AreaResult<RegionContent> areaResult) {
		
		if (areaResult.content.values.get(ColumnType.VALUE) instanceof Integer && 
				((Integer)areaResult.content.values.get(ColumnType.VALUE)) == 128) {
			
			System.out.println("" + areaResult.content.values + areaResult.content.region);
		}
		
		if (areaResult.status.concise == this.isConcised() && 
				areaResult.content.values.get(ColumnType.STRAND) == null || areaResult.content.values.get(ColumnType.STRAND) == getStrand() && 
				areaResult.content.values.get(ColumnType.VALUE) != null &&
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
                getView().getBpRegion().getLength() > minBpLength);
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
