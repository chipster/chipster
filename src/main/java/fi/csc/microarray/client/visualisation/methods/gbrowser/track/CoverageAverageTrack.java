package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;
import java.util.Map.Entry;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.CoverageStorage;

/**
 * Track for showing the coverage of reads. Profile is drawn by calculating
 * the number of nucleotides hitting each basepair location. Should look
 * similar to {@link CoverageEstimateTrack}, but is exact. Also shows where there are
 * large amounts of SNP's as bars chart.
 * 
 * If reverseColor is not null, then strands are visualised separately and
 * SNPs are disabled.
 * 
 */
public class CoverageAverageTrack extends Track { 
	
	private boolean strandSpecificCoverageType;

	private CoverageStorage coverageStorage = new CoverageStorage();
	
	private Collection<Drawable> getAverageDrawables(Strand strand, Color color) {
				
		Collection<Drawable> drawables = getEmptyDrawCollection();

		// Count maximum y coordinate (the bottom of the track)
		int bottomlineY = 0;
	
		int previousValueY = 0;
		int previousEndX = -1;
		
		//Line color is opaque
		Color lineColor = new Color(color.getRGB(), false);
				
		TreeMap<Region, Float> totals = coverageStorage.getTotalAverageCoverage();
				
		for (Entry<Region, Float> entry : totals.entrySet()) {
			
			Region region = entry.getKey();
			
			if (view == null) {
				continue;
			}
			
			float startX = getView().bpToTrackFloat(region.start);
			float endX = getView().bpToTrackFloat(region.end);
			//Round together with position dividends to get the same result than where next block will start
			//int width = (int)(startX + bpWidth) - (int)startX;
			int width = (int)(endX) - (int)startX;
			
			int profileY = 0;
			
			Float coverage = coverageStorage.getAverage(region, strand);					
			
			if (coverage != null) {
				profileY = (int)(float)coverage;
			} else {
				//this totalBase is on the wrong strand
				continue;
			}
			
			int valueY = (int)(bottomlineY + profileY);
			
			drawables.add(new RectDrawable((int)startX, bottomlineY, width,  valueY, color, null));
			
			//Check if there was a gap between profile blocks
			if (previousEndX < (int)startX) {
				
				//End last block with line
				drawables.add(new LineDrawable(previousEndX, bottomlineY, previousEndX,  previousValueY, lineColor));
				
				//Start next line from the bottom
				previousValueY = 0;
			}
			
			//Draw line between height difference of previous and current block
			drawables.add(new LineDrawable((int)startX, previousValueY, (int)startX,  valueY, lineColor));
			//Draw line on top of the current block
			drawables.add(new LineDrawable((int)startX, valueY, (int)startX + width,  valueY, lineColor));									
			
			previousValueY = valueY;
			previousEndX = (int)startX + width;
		}		
		
		//End last block with line
		drawables.add(new LineDrawable(previousEndX, bottomlineY, previousEndX,  previousValueY, lineColor));
		
		return drawables;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();
				
		if (strandSpecificCoverageType) {

			// add drawables of both strands separately
			drawables.addAll(getAverageDrawables(Strand.FORWARD, GBrowserConstants.FORWARD_COLOR));
			drawables.addAll(getAverageDrawables(Strand.REVERSE, GBrowserConstants.REVERSE_COLOR));
			
		} else {
			
			// add drawables according to sum of both strands
			drawables.addAll(getAverageDrawables(Strand.BOTH, GBrowserConstants.getCoverageColor()));
		}				

		return drawables;
	}

	public void processDataResult(DataResult dataResult) {				
		
		coverageStorage.addAverages(dataResult, view.getRequestRegion());		
	}

	@Override
	public int getTrackHeight() {
		return 100;
	}
	
	@Override
	public void defineDataTypes() {
		addDataType(DataType.COVERAGE_AVERAGE);
	}

	/**
	 * @see GBrowserView#drawView
	 */
	@Override
	public boolean canExpandDrawables() {
		return true;
	}

	public void setStrandSpecificCoverageType(boolean b) {
		this.strandSpecificCoverageType = b;
	}
}
