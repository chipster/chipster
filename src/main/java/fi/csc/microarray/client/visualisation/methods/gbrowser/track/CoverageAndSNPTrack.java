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

import fi.csc.microarray.client.visualisation.methods.gbrowser.BaseStorage;
import fi.csc.microarray.client.visualisation.methods.gbrowser.BaseStorage.Base;
import fi.csc.microarray.client.visualisation.methods.gbrowser.BaseStorage.Nucleotide;
import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Track for showing the coverage of reads. Profile is drawn by calculating
 * the number of nucleotides hitting each basepair location. Should look
 * similar to {@link IntensityTrack}, but is exact. Also shows where there are
 * large amounts of SNP's as bars chart.
 * 
 * If reverseColor is not null, then strands are visualised separately and
 * SNPs are disabled.
 * 
 */
public class CoverageAndSNPTrack extends Track { 

	private long maxBpLength;
	private long minBpLength;

	private Color forwardColor;
	private Color reverseColor;

	private boolean highlightSNP = false;

	private BaseStorage theBaseCacheThang = new BaseStorage();
	private DataSource refFile;
	private Collection<RegionContent> refReads = new TreeSet<RegionContent>();
	private ReadpartDataProvider readpartProvider;

	public CoverageAndSNPTrack(View view, DataSource file, ReadpartDataProvider readpartProvider, DataSource refFile, 
			Color forwardColor, Color reverseColor, long minBpLength, long maxBpLength) {
		super(view, file);
		this.forwardColor = forwardColor;
		this.reverseColor = reverseColor;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
		this.readpartProvider = readpartProvider;

		setStrand(Strand.BOTH);
		
		this.refFile = refFile;
	}

	@Override	
	public void initializeListener() {
		super.initializeListener();
		
		// Add listener for reference file
		if (file != null && refFile != null) {
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
		char[] refSeq = SeqBlockTrack.getReferenceArray(refReads, view, Strand.FORWARD);

		// Count nucleotides for each location
		theBaseCacheThang.getNucleotideCounts(readpartProvider.getReadparts(dataStrand), view, refSeq); 

		// Count width of a single bp in pixels
		float bpWidth = ((float) (getView().getWidth()) / getView().getBpRegion().getLength());

		// Count maximum y coordinate (the bottom of the track)
		int bottomlineY = 0;

		// prepare lines that make up the profile for drawing
		Iterator<Base> bases = theBaseCacheThang.iterator();

		// draw lines for each bp region that has some items
		while (bases.hasNext()) {
			Base currentBase = bases.next();

			float startX = getView().bpToTrackFloat(new BpCoord(currentBase.getBpLocation(), chr));
			//Round together with position dividends to get the same result than where next block will start
			int endX = (int)(startX + bpWidth) - (int)startX;
			int profileY = currentBase.getCoverage();
			
			drawables.add(new RectDrawable((int)startX, bottomlineY, endX,  (int)(bottomlineY + profileY), color, null));

			drawSNPBar(drawables, (int)bpWidth, bottomlineY, currentBase, (int)startX);
		}

		return drawables;
	}

	private void drawSNPBar(Collection<Drawable> drawables, int bpWidth, int bottomlineY, Base currentBase, int endX) {
		
		if (highlightSNP && currentBase.hasSignificantSNPs()) {
			int y = bottomlineY;				

			for (Nucleotide nt : Nucleotide.values()) {

				int increment = currentBase.getSNPCounts()[nt.ordinal()];

				if (increment > 0) {
					Color c = SeqBlockTrack.charColors[nt.ordinal()];

					drawables.add(new RectDrawable(endX, y, bpWidth, increment, c, null));

					y += increment;
				}
			}
		}
	}


	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();
		
		if (reverseColor == null) {

			// add drawables according to sum of both strands
			drawables.addAll(getDrawableReads(Strand.BOTH, forwardColor));
			
		} else {
			
			// add drawables of both strands separately
			drawables.addAll(getDrawableReads(Strand.FORWARD, forwardColor));
			drawables.addAll(getDrawableReads(Strand.REVERSED, reverseColor));

		}

		return drawables;
	}

	public void processAreaResult(AreaResult areaResult) {

		// Do not listen to actual read data, because that is taken care by ReadpartDataProvider
		
		// "Spy" on reference sequence data, if available
		if (areaResult.getStatus().file == refFile) {
			this.refReads.addAll(areaResult.getContents());
		}
	}

	@Override
	public Integer getHeight() {
		if (isVisible()) {
			//return super.getHeight();
			
			return 100;
		} else {
			return 0;
		}
	}

	@Override
	public boolean isStretchable() {
		// stretchable unless hidden
		//return isVisible();
		
		return false;
	}

	@Override
	public boolean isVisible() {
		// visible region is not suitable
		return (super.isVisible() &&
				getView().getBpRegion().getLength() > minBpLength &&
				getView().getBpRegion().getLength() <= maxBpLength);
	}

	@Override
	public Map<DataSource, Set<ColumnType>> requestedData() {
		HashMap<DataSource, Set<ColumnType>> datas = new
		HashMap<DataSource, Set<ColumnType>>();
		datas.put(file, new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {
				ColumnType.ID, 
				ColumnType.SEQUENCE,
				ColumnType.STRAND,
				ColumnType.QUALITY,
				ColumnType.CIGAR})));
		
		return datas;
	}

	@Override
	public boolean isConcised() {
		return false;
	}

	/**
	 * @see View#drawView
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
}
