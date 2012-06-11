package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomeBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
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

    final public static int SAMPLING_GRANULARITY = 500;

	private SortedSet<RegionContent> values = new TreeSet<RegionContent>();
	private LinkedList<RegionContent> valueStorageOrder = new LinkedList<RegionContent>();
	private long minBpLength;
	private Color color;
	private boolean doLog;
	private boolean removeTooWide;

	public IntensityTrack(View view, DataSource file, Color c, long maxBpLength, boolean doLog, boolean removeTooWide) {
		super(view, file);
		this.color = c;
		this.doLog = doLog;
		this.minBpLength = maxBpLength;
		this.removeTooWide = removeTooWide;
	}

	@Override
	public Collection<Drawable> getDrawables() {

		Collection<Drawable> drawables = getEmptyDrawCollection();

		Iterator<RegionContent> iterator = values.iterator();
		while (iterator.hasNext()) {

			RegionContent regCont = iterator.next();
			
			// remove values that have gone out of view
			if (!regCont.region.intersects(getView().getBpRegion())) {
				iterator.remove();
				continue;
			}
			
			// remove values that are too wide for this view (when zooming in)
			if (removeTooWide && regCont.region.getLength() > ((getView().getBpRegion().getLength() / SAMPLING_GRANULARITY) * 2)) {
				iterator.remove();
				continue;
			}
			
			// remove values that are too narrow to show (when zooming out)
			if (regCont.region.getLength() < (getView().getBpRegion().getLength() / (SAMPLING_GRANULARITY * 4))) {
				iterator.remove();
				continue;
			}
			
			// do the plotting for this concised value
			int x1 = getView().bpToTrack(regCont.region.start);
			int x2 = getView().bpToTrack(regCont.region.end) + 2;
			int y = 0;						

			double count = (Float) (regCont.values.get(ColumnType.VALUE));
			if (doLog) {
				count = Math.log(count);
			}
			
			int height = (int) Math.min(count * (GenomeBrowserConstants.READ_HEIGHT + GenomeBrowserConstants.SPACE_BETWEEN_READS), getView().getTrackHeight());			

			drawables.add(new RectDrawable(x1, y, x2 - x1, height, color, null));

		}

		// FIXME move this to "gone out of view" place?
		// remove values when they get "too big"
//		while (values.size() > MAX_VALUE_COUNT) {
//			RegionContent oldest = valueStorageOrder.pop();
//			values.remove(oldest);
//		}
		
		return drawables;
	}

	public void processAreaResult(AreaResult areaResult) {		

		for (RegionContent content : areaResult.getContents()) {
			if (areaResult.getStatus().concise == this.isConcised() && 
					content.values.get(ColumnType.STRAND) == getStrand() && 
					content.values.get(ColumnType.VALUE) != null &&
					content.region.intersects(getView().getBpRegion())) { 
				
				values.add(content);
				valueStorageOrder.add(content);
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
		return true;
	}
}
