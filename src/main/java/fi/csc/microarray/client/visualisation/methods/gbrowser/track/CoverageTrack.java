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

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.BaseStorage;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.BaseStorage.Base;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.BaseStorage.Nucleotide;

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

	private long maxBpLength;
	private long minBpLength;

	private boolean highlightSNP = false;

	private BaseStorage theBaseCacheThang = new BaseStorage();
	private AreaRequestHandler refFile;
	private Collection<RegionContent> refReads = new TreeSet<RegionContent>();
	private ReadpartDataProvider readpartProvider;
	private boolean strandSpecificCoverageType;

	public CoverageTrack(ReadpartDataProvider readpartProvider, AreaRequestHandler refFile, long minBpLength, long maxBpLength) {

		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
		this.readpartProvider = readpartProvider;
		
		this.refFile = refFile;
	}

	@Override	
	public void initializeListener() {
		super.initializeListener();
		
		// Add listener for reference file
		if (areaRequestHandler != null && refFile != null) {
			view.getQueueManager().addResultListener(refFile, this);
		}
	}

	/**
	 * Get drawables for a collection of reads.
	 * 
	 * @return
	 */
	private Collection<Drawable> getDrawableReads(Strand dataStrand, Color color) {
		
		Collection<Drawable> drawables = getEmptyDrawCollection();

		Chromosome chr = getView().getBpRegion().start.chr;
		
		// If SNP highlight mode is on, we need reference sequence data
		char[] refSeq = ReadPileTrack.getReferenceArray(refReads, view, Strand.FORWARD);

		// Count nucleotides for each location
		theBaseCacheThang.getNucleotideCounts(readpartProvider.getReadparts(dataStrand), view, refSeq); 

		// Count width of a single bp in pixels
		float bpWidth = (float) (getView().getWidth() / getView().getBpRegionDouble().getLength());

		// Count maximum y coordinate (the bottom of the track)
		int bottomlineY = 0;

		// prepare lines that make up the profile for drawing
		Iterator<Base> bases = theBaseCacheThang.iterator();
		
		int previousValueY = 0;
		int previousEndX = -1;
		
		//Line color is opaque
		Color lineColor = new Color(color.getRGB(), false);

		// draw lines for each bp region that has some items
		while (bases.hasNext()) {
			Base currentBase = bases.next();

			float startX = getView().bpToTrackFloat(new BpCoord(currentBase.getBpLocation(), chr));
			//Round together with position dividends to get the same result than where next block will start
			int width = (int)(startX + bpWidth) - (int)startX;
			int profileY = currentBase.getCoverage();
			
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

			drawSNPBar(drawables, (int)bpWidth, bottomlineY, currentBase, (int)startX);
			
			previousValueY = valueY;
			previousEndX = (int)startX + width;
		}
		
		//End last block with line
		drawables.add(new LineDrawable(previousEndX, bottomlineY, previousEndX,  previousValueY, lineColor));

		return drawables;
	}

	private void drawSNPBar(Collection<Drawable> drawables, int bpWidth, int bottomlineY, Base currentBase, int endX) {
		
		if (highlightSNP && currentBase.hasSignificantSNPs()) {
			int y = bottomlineY;				

			for (Nucleotide nt : Nucleotide.values()) {

				int increment = currentBase.getSNPCounts()[nt.ordinal()];

				if (increment > 0) {
					Color c = GBrowserConstants.charColors[nt.ordinal()];

					drawables.add(new RectDrawable(endX, y, bpWidth, increment, c, null));

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
			drawables.addAll(getDrawableReads(Strand.FORWARD, GBrowserConstants.FORWARD_COLOR));
			drawables.addAll(getDrawableReads(Strand.REVERSE, GBrowserConstants.REVERSE_COLOR));
			
		} else {
			
			// add drawables according to sum of both strands
			drawables.addAll(getDrawableReads(Strand.BOTH, GBrowserConstants.COVERAGE_COLOR));
		}

		return drawables;
	}

	public void processAreaResult(AreaResult areaResult) {

		// Do not listen to actual read data, because that is taken care by ReadpartDataProvider
		
		// "Spy" on reference sequence data, if available
		if (areaResult.getStatus().areaRequestHandler == refFile) {
			this.refReads.addAll(areaResult.getContents());
		}
	}

	@Override
	public int getHeight() {
		return 100;
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
		datas.put(areaRequestHandler, new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {ColumnType.COVERAGE}))); 
		
		// We might also need reference sequence data
		if (highlightSNP && this.getView().getBpRegion().getLength() < this.getView().getWidth() * 2) {
			datas.put(refFile, new HashSet<ColumnType>(Arrays.asList(new ColumnType[] { ColumnType.SEQUENCE })));
		}
		
		return datas;
	}

	/**
	 * @see GBrowserView#drawView
	 */
	@Override
	public boolean canExpandDrawables() {
		return true;
	}

	public void enableSNPHighlight() {
		// turn on highlighting mode
		highlightSNP = true;
	}

	public void disableSNPHighlight() {
		highlightSNP = false;
	}

	public void setStrandSpecificCoverageType(boolean b) {
		this.strandSpecificCoverageType = b;
	}
}
