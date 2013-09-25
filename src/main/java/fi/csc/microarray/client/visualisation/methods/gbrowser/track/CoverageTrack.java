package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;
import java.util.TreeMap;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.BaseStorage.Base;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.BaseStorage.Nucleotide;
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
public class CoverageTrack extends Track { 

	private boolean highlightSNP = false;

	private Collection<Feature> refReads = new TreeSet<Feature>();
	private boolean strandSpecificCoverageType;
	private Integer detailsIndex = null;
	private Integer referenceIndex = null;

	private CoverageStorage coverageStorage = new CoverageStorage();

	public CoverageTrack(DataThread coverage, DataThread referenceSequenceFile) {
		
		super();

		detailsIndex = addDataThread(coverage);

		if (referenceSequenceFile != null) {
			referenceIndex = addDataThread(referenceSequenceFile);
		}
	}

	@Override	
	public void initializeListener() {
		super.initializeListener();
	}

	private Collection<Drawable> getCoverageDrawables(Strand strand, Color color) {				
		
		Collection<Drawable> drawables = getEmptyDrawCollection();

		Chromosome chr = getView().getBpRegion().start.chr;
		
		// If SNP highlight mode is on, we need reference sequence data
		char[] refSeq = ReadPileTrack.getReferenceArray(refReads, view, Strand.FORWARD);

		// Count width of a single bp in pixels
		float bpWidth = (float) (getView().getWidth() / getView().getBpRegionDouble().getLength());

		// Count maximum y coordinate (the bottom of the track)
		int bottomlineY = 0;
	
		int previousValueY = 0;
		int previousEndX = -1;
		
		//Line color is opaque
		Color lineColor = new Color(color.getRGB(), false);
				
		TreeMap<BpCoord, Base> totalBases = coverageStorage.getTotalBases();
		
		for (Base base : totalBases.values()) {
			
			Nucleotide reference = null;
			int viewIndex = (int) (base.getBpLocation() - view.getBpRegion().start.bp);
			if (viewIndex >= 0 && viewIndex < refSeq.length) {
				reference = Nucleotide.fromCharacter(refSeq[viewIndex]);
				Base refBase = new Base(base.getBpLocation(), reference);
				refBase.setNucleotideCounts(base.getNucleotideCounts());
				base = refBase;
			}
						
			BpCoord location = new BpCoord(base.getBpLocation(), chr);
			float startX = getView().bpToTrackFloat(location);
			//Round together with position dividends to get the same result than where next block will start
			int width = (int)(startX + bpWidth) - (int)startX;
			
			int profileY = 0;
			
			Base coverageBase = coverageStorage.getBase(location, strand);			
			
			if (coverageBase != null) {
				profileY = coverageBase.getCoverage();
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
			
			drawSNPBar(drawables, (int)bpWidth, bottomlineY, base, strand, (int)startX);
			
			previousValueY = valueY;
			previousEndX = (int)startX + width;
		}		
		
		//End last block with line
		drawables.add(new LineDrawable(previousEndX, bottomlineY, previousEndX,  previousValueY, lineColor));

		return drawables;
	}

	private void drawSNPBar(Collection<Drawable> drawables, int bpWidth, int bottomlineY, Base base, Strand strand, int endX) {
		
		if (!strandSpecificCoverageType && highlightSNP && base.hasSignificantSNPs()) {
			
			int y = bottomlineY;				

			for (Nucleotide nt : Nucleotide.values()) {
				
				int increment = 0;				
				increment += base.getSNPCounts()[nt.ordinal()];

				if (increment > 0) {
					Color c = GBrowserConstants.charColors[nt.ordinal()];

					drawables.add(new RectDrawable(endX, y, bpWidth, increment, c, Color.black));

					y += increment;
				}
			}
		}
	}


	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();
				
		if (strandSpecificCoverageType) {

			// add drawables of both strands separately
			drawables.addAll(getCoverageDrawables(Strand.FORWARD, GBrowserConstants.FORWARD_COLOR));
			drawables.addAll(getCoverageDrawables(Strand.REVERSE, GBrowserConstants.REVERSE_COLOR));
			
		} else {
			
			// add drawables according to sum of both strands
			drawables.addAll(getCoverageDrawables(Strand.BOTH, GBrowserConstants.getCoverageColor()));
		}				

		return drawables;
	}

	public void processDataResult(DataResult dataResult) {				
		
		coverageStorage.addBaseCoverage(dataResult, view.getRequestRegion());
		
		// "Spy" on reference sequence data, if available
		if (dataResult.getStatus().getDataThread() == dataThreads.get(referenceIndex)) {
			this.refReads.addAll(dataResult.getFeatures());
		}
	}

	@Override
	public int getTrackHeight() {
		return 100;
	}
	
	@Override
	public void defineDataTypes() {
		addDataType(dataThreads.get(detailsIndex), DataType.COVERAGE);
		
		// We might also need reference sequence data
		if (highlightSNP && this.getView().getBpRegion().getLength() < this.getView().getWidth() * 2) {
			addDataType(dataThreads.get(referenceIndex), DataType.SEQUENCE);
		}
	}

	/**
	 * @see GBrowserView#drawView
	 */
	@Override
	public boolean canExpandDrawables() {
		return true;
	}

	public void setSNPHighlight(boolean highlightSnp) {
		this.highlightSNP = highlightSnp;
	}

	public void setStrandSpecificCoverageType(boolean b) {
		this.strandSpecificCoverageType = b;
	}
}
