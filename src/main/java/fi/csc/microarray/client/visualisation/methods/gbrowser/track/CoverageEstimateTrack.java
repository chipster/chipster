package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedSet;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Generic track for showing high level distribution of items (genes, transcripts, reads...) on the genome.
 * The result is an approximation.
 *
 */
public class CoverageEstimateTrack extends Track {

    final public static int SAMPLING_GRANULARITY = 10;

	private static final int MAX_VALUE_COUNT = 1000;

	private SortedSet<RegionContent> values = new TreeSet<RegionContent>();
	private LinkedList<RegionContent> valueStorageOrder = new LinkedList<RegionContent>();
	private long minBpLength;

	private boolean strandSpecificCoverageType;

	public CoverageEstimateTrack(long maxBpLength) {
		
		this.minBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {

		Collection<Drawable> drawables = getEmptyDrawCollection();
	
//		Color bg = new Color(0f, 0f, 0f, 0.05f);
//		drawables.add(new RectDrawable(0,  0,  getView().getWidth(), getHeight(), bg, bg));
		
		// remove values when they get "too big"
		while (values.size() > MAX_VALUE_COUNT) {
			RegionContent oldest = valueStorageOrder.pop();
			values.remove(oldest);
		}
		
		List<RegionValue> forward = new LinkedList<RegionValue>();
		List<RegionValue> reverse = new LinkedList<RegionValue>();		
		
		Iterator<RegionContent> iterator = values.iterator();
		while (iterator.hasNext()) {

			RegionContent regCont = iterator.next();
			
			// remove values that have gone out of view
			if (!getView().requestIntersects(regCont.region)) {
				
				iterator.remove();
				continue;
			}
						
			int x1 = getView().bpToTrack(regCont.region.start);
			int x2 = getView().bpToTrack(regCont.region.end);
			
			x2 = Math.max(x2, x1 + 2);
			
			int fCount = (Integer) (regCont.values.get(ColumnType.COVERAGE_ESTIMATE_FORWARD));
			int rCount = (Integer) (regCont.values.get(ColumnType.COVERAGE_ESTIMATE_REVERSE));							
			
			if (!strandSpecificCoverageType) {				
				fCount += rCount;
				rCount = 0;
			}

			forward.add(new RegionValue(x1, x2, fCount / (float)regCont.region.getLength()));
			reverse.add(new RegionValue(x1, x2, rCount / (float)regCont.region.getLength()));
		}
			
		float[] continuousForward = makeContinuous(forward, getView().getWidth(), getView().getWidth() / SAMPLING_GRANULARITY * 4);
		float[] continuousReverse = makeContinuous(reverse, getView().getWidth(), getView().getWidth() / SAMPLING_GRANULARITY * 4);		
		
		int[] smoohtForward = smooth(continuousForward);
		int[] smoohtReverse = smooth(continuousReverse);
		
		int y = 0;
		
		Color forwardColor;
		Color reverseColor;
		
		if (strandSpecificCoverageType) {
			forwardColor = GBrowserConstants.FORWARD_COLOR;
			reverseColor = GBrowserConstants.REVERSE_COLOR;
		} else {
			forwardColor = GBrowserConstants.COVERAGE_COLOR;
			reverseColor = null;
		}
			
		for (int i = 0; i < smoohtForward.length; i++) {
			
			int fValue = (int) smoohtForward[i];			
			int rValue = (int) smoohtReverse[i];
									
			drawables.add(new RectDrawable(i, y, 1, fValue, forwardColor, null));
			
			if (strandSpecificCoverageType) {
				drawables.add(new RectDrawable(i, y, 1, rValue, reverseColor, null));
			}
		}
		
		return drawables;
	}
	
	private class RegionValue {
		public RegionValue(int start, int end, float value) {
			this.start = start;
			this.end = end;
			this.value = value;
		}
		int start;
		int end;
		float value;
	}

	private float[] makeContinuous(List<RegionValue> points, int width, int maxPointDistance) {
		
		RegionValue lastPoint = null;
		
		float[] continuous = new float[width];
		
		for (RegionValue point : points) {
					
			for (int i = point.start; i < point.end; i++) {
				if (i >= 0 && i < width) {
					continuous[i] = point.value;  
				}
			}				

			if (lastPoint != null) {

					if (lastPoint.end < point.start && point.start - lastPoint.start <= maxPointDistance) {
	
						for (int i = lastPoint.end; i < point.start; i++) {
							if (i >= 0 && i < width) {
								continuous[i] = (point.value + lastPoint.value) / 2;  
							}
						}
					}
			}
			
			lastPoint = point;
		}
		
		return continuous;		
	}

	private int[] smooth(float[] values) {
		
		final int WINDOW = 8;
		
		int[] smooth = new int[values.length];
		
		for (int i = 0; i < values.length; i++) {
			int sum = 0;
			int divisor = WINDOW;
			for (int j = i - WINDOW / 2; j < i + WINDOW / 2; j++) {
				if (j >= 0 && j < values.length) {
					sum += values[j];
				} else {
					divisor--;
				}
			}
			smooth[i] = sum / divisor;
			//smooth[i] = (int) values[i]; //turn off smoothing
		}
				
		return smooth;
	}

	public void processAreaResult(AreaResult areaResult) {		

		for (RegionContent content : areaResult.getContents()) {
			if (getView().requestIntersects(content.region) && content.values.containsKey(ColumnType.COVERAGE_ESTIMATE_FORWARD)) {
								
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
	
	/**
	 * @see GBrowserView#drawView
	 */
	@Override
	public boolean canExpandDrawables() {
		return true;
	}

	public void setStrandSpecificCoverageType(boolean b) {
		strandSpecificCoverageType = b;
	}
}
