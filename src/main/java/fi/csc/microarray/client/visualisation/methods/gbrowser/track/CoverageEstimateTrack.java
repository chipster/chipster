package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;
import java.util.SortedSet;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;

/**
 * Generic track for showing high level distribution of items (genes, transcripts, reads...) on the genome.
 * The result is an approximation.
 *
 */
public class CoverageEstimateTrack extends ScaleTrack {

    public CoverageEstimateTrack() {
		super(100, 3);
	}

	final public static int SAMPLING_GRANULARITY = 4;

	private static final int MAX_VALUE_COUNT = 1000;

	private SortedSet<Feature> values = new TreeSet<Feature>();
	private LinkedList<Feature> valueStorageOrder = new LinkedList<Feature>();

	private boolean strandSpecificCoverageType;

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();
		
		super.updateScale();
		
		// don't show scale because these are only estimates (but it is still used for scaling)
		//drawables.addAll(getScaleDrawables());		
		drawables.addAll(getEstimateDrawables());
							
		return drawables;
	}
		
	public int getMaxTotalCoverage() {
		float max = 0;
		
		for (Feature feature : values) {
			
			float forward = ((Integer) (feature.values.get(DataType.COVERAGE_ESTIMATE_FORWARD)) / (float)feature.region.getLength());
			float reverse = ((Integer) (feature.values.get(DataType.COVERAGE_ESTIMATE_REVERSE)) / (float)feature.region.getLength());
			
			float estimate = forward + reverse;
			
			if (estimate > max){
				max = estimate;
			}
		}
		return (int) max;
	}

	public Collection<Drawable> getEstimateDrawables() {

		Collection<Drawable> drawables = getEmptyDrawCollection();
		
		// remove values when they get "too big"
		while (values.size() > MAX_VALUE_COUNT) {
			Feature oldest = valueStorageOrder.pop();
			values.remove(oldest);
		}
		
		List<RegionValue> forward = new LinkedList<RegionValue>();
		List<RegionValue> reverse = new LinkedList<RegionValue>();		
		
		Iterator<Feature> iterator = values.iterator();
		while (iterator.hasNext()) {

			Feature regCont = iterator.next();
			
			// remove values that have gone out of view
			if (!getView().requestIntersects(regCont.region)) {
				
				iterator.remove();
				continue;
			}
						
			int x1 = getView().bpToTrack(regCont.region.start);
			int x2 = getView().bpToTrack(regCont.region.end);
			
			x2 = Math.max(x2, x1 + 2);
			
			int fCount = (Integer) (regCont.values.get(DataType.COVERAGE_ESTIMATE_FORWARD));
			int rCount = (Integer) (regCont.values.get(DataType.COVERAGE_ESTIMATE_REVERSE));							
			
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
			forwardColor = GBrowserConstants.getCoverageColor();
			reverseColor = null;
		}
			
		for (int i = 0; i < smoohtForward.length; i++) {
			
			int fValue = super.getScaledY(smoohtForward[i]);			
			int rValue = super.getScaledY(smoohtReverse[i]);
									
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

	public void processDataResult(DataResult dataResult) {		

		for (Feature content : dataResult.getFeatures()) {
			if (getView().requestIntersects(content.region) && content.values.containsKey(DataType.COVERAGE_ESTIMATE_FORWARD)) {
								
				values.add(content);
				valueStorageOrder.add(content);
			}
		}
	}   
    
    @Override
	public void defineDataTypes() {
		
		if (isSuitableViewLength()) {
			
			addDataType(DataType.COVERAGE_ESTIMATE_FORWARD);
			addDataType(DataType.COVERAGE_ESTIMATE_REVERSE);
			
		} else {
			
			addDataType(DataType.CANCEL);
		}
	}

	public void setStrandSpecificCoverageType(boolean b) {
		strandSpecificCoverageType = b;
	}
}
