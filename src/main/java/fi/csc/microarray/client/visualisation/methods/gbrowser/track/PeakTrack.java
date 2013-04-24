package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Track for showing the location of predicted peaks. Peaks cannot overlap. 
 *
 */
public class PeakTrack extends Track {

	private static final int MIN_VISIBLE_PEAK_SIZE = 5;

	private static final int PEAK_SYMBOL_HEIGHT = 10;

	private Collection<RegionContent> peaks = new TreeSet<RegionContent>();

	private long maxBpLength;
	private long minBpLength;

	private Color color;


	public PeakTrack(Color color, long minBpLength, long maxBpLength) {

		this.color = color;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		if (peaks != null) {

			Iterator<RegionContent> iter = peaks.iterator();
			while (iter.hasNext()) {

				RegionContent peak = iter.next();

				if (!getView().requestIntersects(peak.region)) {
					iter.remove();
					continue;
				}

				createDrawable(peak.region.start, peak.region.end, PEAK_SYMBOL_HEIGHT, color, drawables);
			}
		}

		return drawables;
	}

	private void createDrawable(BpCoord startBp, BpCoord endBp, int height, Color c, Collection<Drawable> drawables) {
		Rectangle rect = new Rectangle();

		rect.x = getView().bpToTrack(startBp);
		rect.width = getView().bpToTrack(endBp) - rect.x;
		
		if (rect.width < MIN_VISIBLE_PEAK_SIZE) {
			rect.width = MIN_VISIBLE_PEAK_SIZE;
		}

		rect.y = getHeight() / 2;
		rect.height = height;

		drawables.add(new RectDrawable(rect, c, c.darker()));
	}

	public void processAreaResult(AreaResult areaResult) {

		this.peaks.addAll(areaResult.getContents());
		getView().redraw();
	}
    
    @Override
    public boolean isVisible() {
        // visible region is not suitable
        return (super.isVisible() &&
                getView().getBpRegion().getLength() > minBpLength &&
                getView().getBpRegion().getLength() <= maxBpLength);
    }
	
    @Override
    public Map<AreaRequestHandler, Set<ColumnType>> requestedData() {
        HashMap<AreaRequestHandler, Set<ColumnType>> datas = new
        HashMap<AreaRequestHandler, Set<ColumnType>>();
        datas.put(areaRequestHandlers.get(0), new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {
                ColumnType.CHROMOSOME,
                ColumnType.BP_START,
                ColumnType.BP_END })));
        return datas;
    }
	
	@Override
	public int getHeight() {
		return PEAK_SYMBOL_HEIGHT * 2;
	}
}
