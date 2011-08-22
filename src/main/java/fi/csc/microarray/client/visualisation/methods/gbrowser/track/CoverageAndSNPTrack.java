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

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.BaseStorage;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.BaseStorage.Nucleotide;
import fi.csc.microarray.client.visualisation.methods.gbrowser.BaseStorage.Base;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
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
 * similar to IntensityTrack, but is exact. Also shows where there are
 * large amounts of SNP's as bars chart.
 * 
 *
 * @see IntensityTrack
 */
public class CoverageAndSNPTrack extends Track {

	private long maxBpLength;
	private long minBpLength;

	private Color forwardColor;

	private boolean highlightSNP = false;

	private BaseStorage theBaseCacheThang = new BaseStorage();
	private DataSource refFile;
	private Collection<RegionContent> refReads = new TreeSet<RegionContent>();
	private ReadpartDataProvider readpartProvider;

	public CoverageAndSNPTrack(View view, DataSource file, ReadpartDataProvider readpartProvider, Class<? extends AreaRequestHandler> handler, DataSource refFile, Class<? extends AreaRequestHandler> refHandler, 
			Color forwardColor, long minBpLength, long maxBpLength) {
		super(view, file, handler);
		this.forwardColor = forwardColor;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
		this.readpartProvider = readpartProvider;

		setStrand(Strand.BOTH);
		
		this.refFile = refFile;
		view.getQueueManager().createQueue(refFile, refHandler);
	}

	@Override	
	public void initializeListener() {
		super.initializeListener();
		
		// Add listener for reference file
		if (file != null) {
			view.getQueueManager().addResultListener(refFile, this);
		}
	}

	/**
	 * Get drawables for a collection of reads.
	 * 
	 * @return
	 */
	private Collection<Drawable> getDrawableReads(Color color) {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		Chromosome chr = getView().getBpRegion().start.chr;

		// If SNP highlight mode is on, we need reference sequence data
		char[] refSeq = SeqBlockTrack.getReferenceArray(refReads, view, strand);

		// Count nucleotides for each location
		theBaseCacheThang.getNucleotideCounts(readpartProvider.getReadparts(getStrand()), view, refSeq); 

		// Count width of a single bp in pixels
		int bpWidth = (int) (getView().getWidth() / getView().getBpRegion().getLength());

		// Count maximum y coordinate (the bottom of the track)
		int bottomlineY = 0;

		// prepare lines that make up the profile for drawing
		Iterator<Base> bases = theBaseCacheThang.iterator();
		if (bases.hasNext()) {

			Base previousBase = bases.next();

			// draw a line from the beginning of the graph to the first location
			int startX = getView().bpToTrack(new BpCoord(previousBase.getBpLocation(), chr));
			long startY = previousBase.getCoverage();
			drawables.add(new LineDrawable(0, bottomlineY,
					(int)(startX - bpWidth), bottomlineY, color));
			drawables.add(new LineDrawable((int)(startX - bpWidth), bottomlineY,
					startX, (int)(bottomlineY + startY), color));

			// Draw bar for the first base
			drawSNPBar(drawables, bpWidth, bottomlineY, previousBase, startX);

			// draw lines for each bp region that has some items
			while (bases.hasNext()) {
				Base currentBase = bases.next();
				
				startX = getView().bpToTrack(new BpCoord(previousBase.getBpLocation(), chr));
				startY = previousBase.getCoverage();
				int endX = getView().bpToTrack(new BpCoord(currentBase.getBpLocation(), chr));
				long endY = currentBase.getCoverage();

				drawSNPBar(drawables, bpWidth, bottomlineY, currentBase, endX);

				if (currentBase.getBpLocation() - previousBase.getBpLocation() == 1) {
					// join adjacent bp locations with a line
					drawables.add(new LineDrawable(startX, (int)(bottomlineY + startY),
							endX, (int)(bottomlineY + endY), color));
				} else {
					// Join locations that are more than one bp apart using 3 lines
					// 1. Slope down
					drawables.add(new LineDrawable((int)startX, (int)(bottomlineY + startY),
							(int)(startX + bpWidth), bottomlineY, color));
					// 2. Bottomline level for the empty area
					drawables.add(new LineDrawable((int)(startX + bpWidth), bottomlineY,
							(int)(endX - bpWidth), bottomlineY, color));
					// 3. Slope up
					drawables.add(new LineDrawable((int)(endX - bpWidth), bottomlineY,
							(int)endX, (int)(bottomlineY + endY), color));
				}
				

				previousBase = currentBase;
			}

			// Draw a line from the last location to the end of the graph
			int endX = getView().bpToTrack(new BpCoord(previousBase.getBpLocation(), chr));
			long endY = previousBase.getCoverage();
			drawables.add(new LineDrawable(endX, (int)(bottomlineY + endY),
					(int)(endX + bpWidth), bottomlineY, color));
			drawables.add(new LineDrawable((int)(endX + bpWidth), bottomlineY,
					getView().getWidth(), bottomlineY, color));
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

					drawables.add(new RectDrawable(endX, y, bpWidth, increment, c, c));

					y += increment;
				}
			}
		}
	}


	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		// add drawables from both reads (if present)
		drawables.addAll(getDrawableReads(forwardColor));

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

	@Override
	public String getName() {
		return "ProfileSNPTrack";
	}

	public void enableSNPHighlight() {
		// turn on highlighting mode
		highlightSNP = true;
	}

	public void disableSNPHighlight() {
		highlightSNP = false;
	}
}
