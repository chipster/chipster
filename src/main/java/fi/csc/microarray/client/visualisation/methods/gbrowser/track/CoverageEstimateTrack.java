package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Generic track for showing high level distribution of items (genes, transcripts, reads...) on the genome.
 * The result is an approximation.
 *
 */
public class CoverageEstimateTrack extends Track {

    final public static int SAMPLING_GRANULARITY = 100;

	private static final int MAX_VALUE_COUNT = SAMPLING_GRANULARITY * 2;

	private SortedSet<RegionContent> values = new TreeSet<RegionContent>();
	private LinkedList<RegionContent> valueStorageOrder = new LinkedList<RegionContent>();
	private long minBpLength;
	private Color color;
	private boolean doLog;
	private boolean removeTooWide;

	public CoverageEstimateTrack(Color c, long maxBpLength, boolean doLog, boolean removeTooWide) {

		this.color = c;
		this.doLog = doLog;
		this.minBpLength = maxBpLength;
		this.removeTooWide = removeTooWide;
	}

	@Override
	public Collection<Drawable> getDrawables() {

		Collection<Drawable> drawables = getEmptyDrawCollection();
		
		
		// remove values when they get "too big"
		while (values.size() > MAX_VALUE_COUNT) {
			RegionContent oldest = valueStorageOrder.pop();
			values.remove(oldest);
		}

		Iterator<RegionContent> iterator = values.iterator();
		while (iterator.hasNext()) {

			RegionContent regCont = iterator.next();
			
			// remove values that have gone out of view
			if (!getView().requestIntersects(regCont.region)) {
				iterator.remove();
				continue;
			}
			
			// remove values that are too wide for this view (when zooming in)
//			if (removeTooWide && regCont.region.getLength() > ((getView().getBpRegion().getLength() / SAMPLING_GRANULARITY) * 2)) {
//				iterator.remove();
//				continue;
//			}
			
//			// remove values that are too narrow to show (when zooming out)
//			if (regCont.region.getLength() < (getView().getBpRegion().getLength() / (SAMPLING_GRANULARITY * 4))) {
//				iterator.remove();
//				continue;
//			}
			
			// do the plotting for this concised value
			int x1 = getView().bpToTrack(regCont.region.start);
			int x2 = getView().bpToTrack(regCont.region.end) + 2;
			int y = 0;						
			
			double count;
			if (this.getStrand() == Strand.FORWARD) {				
				count = (Integer) (regCont.values.get(ColumnType.COVERAGE_ESTIMATE_FORWARD));
			} else  {
				count = (Integer) (regCont.values.get(ColumnType.COVERAGE_ESTIMATE_REVERSE));
			}
			
			if (doLog) {
				count = Math.log(count);
			}
			
			int height = (int) Math.min(count * (2), getHeight());			

			drawables.add(new RectDrawable(x1, y, x2 - x1, height, color, null));

		}
		
		return drawables;
	}

	public void processAreaResult(AreaResult areaResult) {		

		for (RegionContent content : areaResult.getContents()) {
			if (content.region.intersects(getView().getBpRegion()) && content.values.containsKey(ColumnType.COVERAGE_ESTIMATE_FORWARD)) {
								
				values.add(content);
				valueStorageOrder.add(content);
			}
		}

		getView().redraw();
	}
    
    @Override
    public boolean isVisible() {
        // visible region is not suitable
        return (super.isVisible() &&
                getView().getBpRegion().getLength() > minBpLength);
    }

    @Override
    public Map<AreaRequestHandler, Set<ColumnType>> requestedData() {
        HashMap<AreaRequestHandler, Set<ColumnType>> datas = new
                HashMap<AreaRequestHandler, Set<ColumnType>>();
        datas.put(areaRequestHandler, new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {
				ColumnType.COVERAGE_ESTIMATE_FORWARD,
				ColumnType.COVERAGE_ESTIMATE_REVERSE})));
        return datas;
    }
	
	@Override
	public int getHeight() {
	    return 100;
	}
}
