package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomeBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Cigar;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ReadPart;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.utils.Sequence;

/**
 * Track that shows actual content of reads using color coding.
 * 
 */
public class SeqBlockTrack extends Track {

	public static final Color[] charColors = new Color[] { 
		new Color(64, 192, 64, 128), // A
		new Color(64, 64, 192, 128), // C
		new Color(128, 128, 128, 128), // G
		new Color(192, 64, 64, 128) // T
	};

	private Collection<RegionContent> reads = new TreeSet<RegionContent>();
	private List<Integer> occupiedSpace = new ArrayList<Integer>();

	private long maxBpLength;
	private long minBpLength;

	private DataSource refData;
	private Collection<RegionContent> refReads = new TreeSet<RegionContent>();

	private boolean highlightSNP = false;

	public SeqBlockTrack(View view, DataSource file, Class<? extends AreaRequestHandler> handler, Color fontColor, long minBpLength, long maxBpLength) {
		super(view, file, handler);
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		occupiedSpace.clear();

		// If SNP highlight mode is on, we need reference sequence data
		char[] refSeq = highlightSNP ? getReferenceArray(refReads, view, strand) : null;

		// Main loop: Iterate over RegionContent objects (one object corresponds to one read)
		Iterator<RegionContent> iter = reads.iterator();
		while (iter.hasNext()) {

			RegionContent read = iter.next();

			// Remove reads that are not in this view
			if (!read.region.intercepts(getView().getBpRegion())) {
				iter.remove();
				continue;
			}

			// Collect relevant data for this read

			// Split read into continuous blocks (elements) by using the cigar
			List<ReadPart> visibleRegions = Cigar.splitVisibleElements(read);
			for (ReadPart visibleRegion : visibleRegions) {

				// Skip elements that are not in this view
				if (!visibleRegion.intercepts(getView().getBpRegion())) {
					continue;
				}
				
				// Width in basepairs
				long widthInBps = visibleRegion.getLength();

				// Create rectangle covering the correct screen area (x-axis)
				Rectangle rect = new Rectangle();
				rect.x = getView().bpToTrack(visibleRegion.start);
				rect.width = (int) Math.round(getView().bpWidth() * widthInBps);

				// Do not draw invisible rectangles
				if (rect.width < 2) {
					rect.width = 2;
				}

				// Read parts are drawn in order and placed in layers
				int layer = 0;
				while (occupiedSpace.size() > layer && occupiedSpace.get(layer) > rect.x + 1) {
					layer++;
				}

				// Read part reserves the space of the layer from end to left corner of the screen
				int end = rect.x + rect.width;
				if (occupiedSpace.size() > layer) {
					occupiedSpace.set(layer, end);
				} else {
					occupiedSpace.add(end);
				}

				// Now we can decide the y coordinate
				rect.y = getYCoord(layer, GenomeBrowserConstants.READ_HEIGHT);
				rect.height = GenomeBrowserConstants.READ_HEIGHT;

				// Check if we are about to go over the edge of the drawing area
				boolean lastBeforeMaxStackingDepthCut = getYCoord(layer + 1, GenomeBrowserConstants.READ_HEIGHT) < 0;

				// Check if we are over the edge of the drawing area
				if (rect.y < 0) {
					continue;
				}

				// Check if we have enough space for the actual sequence (at least pixel per nucleotide)
				String seq = visibleRegion.getSequencePart();
				if (rect.width < seq.length()) {
					// Too little space - only show one rectangle for each read part

					Color color = Color.gray;
					
					// Mark last line that will be drawn
					if (lastBeforeMaxStackingDepthCut) {
						color = color.brighter();
					}

					drawables.add(new RectDrawable(rect, color, null));

				} else {
					// Enough space - show color coding for each nucleotide

					// Complement the read if on reverse strand
					if ((Strand) read.values.get(ColumnType.STRAND) == Strand.REVERSED) {

						StringBuffer buf = new StringBuffer(seq.toUpperCase());

						// Complement
						seq = buf.toString().replace('A', 'x'). // switch A and T
						replace('T', 'A').replace('x', 'T').

						replace('C', 'x'). // switch C and G
						replace('G', 'C').replace('x', 'G');
					}

					// Prepare to draw single nucleotides
					float increment = getView().bpWidth();
					float startX = getView().bpToTrackFloat(visibleRegion.start);

					// Draw each nucleotide
					for (int j = 0; j < seq.length(); j++) {

						char letter = seq.charAt(j);

						long refIndex = j;

						// Choose a color depending on viewing mode
						Color bg = Color.white;
						long posInRef = read.region.start.bp.intValue() + refIndex - getView().getBpRegion().start.bp.intValue();
						if (highlightSNP && posInRef >= 0 && posInRef < refSeq.length && Character.toLowerCase(refSeq[(int)posInRef]) == Character.toLowerCase(letter)) {
							bg = Color.gray;
						} else {
							switch (letter) {
							case 'A':
								bg = charColors[0];
								break;
							case 'C':
								bg = charColors[1];
								break;
							case 'G':
								bg = charColors[2];
								break;
							case 'T':
								bg = charColors[3];
								break;
							}
						}

						// Tell that we have reached max. stacking depth
						if (lastBeforeMaxStackingDepthCut) {
							bg = bg.brighter();
						}

						// Draw rectangle
						int x = Math.round(startX + Math.round(((float)refIndex) * increment));
						int width = increment >= 1.0f ? Math.round(increment) : 1;  
						drawables.add(new RectDrawable(x, rect.y, width, GenomeBrowserConstants.READ_HEIGHT, bg, null));
					}
				}
			}
		}

		return drawables;
	}

	private int getYCoord(int layer, int height) {
		return (int) (getView().getTrackHeight() - ((layer + 1) * (height + GenomeBrowserConstants.SPACE_BETWEEN_READS)));
	}

	public void processAreaResult(AreaResult<RegionContent> areaResult) {

		// Check that areaResult has same concised status (currently always false) and correct strand
		if (areaResult.status.file == file && areaResult.status.concise == isConcised() && areaResult.content.values.get(ColumnType.STRAND) == getStrand()) {
			// Add this to queue of RegionContents to be processed
			this.reads.add(areaResult.content);
			getView().redraw();
		}

		// "Spy" on reference sequence data, if available
		if (areaResult.status.file == refData) {
			this.refReads.add(areaResult.content);
		}
	}

	@Override
	public Integer getHeight() {
		if (isVisible()) {
			return super.getHeight();
		} else {
			return 0;
		}
	}

	@Override
	public boolean isStretchable() {
		return isVisible(); // Stretchable unless hidden
	}

	@Override
	public boolean isVisible() {
		// visible region is not suitable
		return (super.isVisible() && getView().getBpRegion().getLength() > minBpLength && getView().getBpRegion().getLength() <= maxBpLength);
	}

	@Override
	public Map<DataSource, Set<ColumnType>> requestedData() {
		HashMap<DataSource, Set<ColumnType>> datas = new HashMap<DataSource, Set<ColumnType>>();
		datas.put(file, new HashSet<ColumnType>(Arrays.asList(new ColumnType[] { ColumnType.SEQUENCE, ColumnType.STRAND, ColumnType.CIGAR })));

		// We might also need reference sequence data
		if (highlightSNP) {
			datas.put(refData, new HashSet<ColumnType>(Arrays.asList(new ColumnType[] { ColumnType.SEQUENCE })));
		}

		return datas;
	}

	@Override
	public boolean isConcised() {
		return false;
	}

	/**
	 * Enable SNP highlighting and set reference data.
	 * 
	 * @param highlightSNP
	 * @see SeqBlockTrack.setReferenceSeq
	 */
	public void enableSNPHighlight(DataSource file, Class<? extends AreaRequestHandler> handler) {
		// turn on highlighting mode
		highlightSNP = true;

		// set reference data
		refData = file;
		view.getQueueManager().createQueue(file, handler);
		view.getQueueManager().addResultListener(file, this);
	}

	/**
	 * Disable SNP highlighting.
	 * 
	 * @param file
	 */
	public void disableSNPHiglight(DataSource file) {
		// turn off highlighting mode
		highlightSNP = false;
	}

	/**
	 * Convert reference sequence reads to a char array.
	 */
	public static char[] getReferenceArray(Collection<RegionContent> refReads, View view, Strand strand) {
		char[] refSeq = new char[0];
		Iterator<RegionContent> iter = refReads.iterator();
		refSeq = new char[view.getBpRegion().getLength().intValue() + 1];
		int startBp = view.getBpRegion().start.bp.intValue();
		int endBp = view.getBpRegion().end.bp.intValue();
		RegionContent read;
		while (iter.hasNext()) {
			read = iter.next();
			if (!read.region.intercepts(view.getBpRegion())) {
				iter.remove();
				continue;
			}

			// we might need to reverse reference sequence
			char[] readBases;
			if (strand == Strand.REVERSED) {
				readBases = Sequence.complement((String) read.values.get(ColumnType.SEQUENCE)).toCharArray();
			} else {
				readBases = ((String) read.values.get(ColumnType.SEQUENCE)).toCharArray();
			}

			int readStart = read.region.start.bp.intValue();
			int readNum = 0;
			int nextPos = 0;
			for (char c : readBases) {
				nextPos = readStart + readNum++;
				if (nextPos >= startBp && nextPos <= endBp) {
					refSeq[nextPos - startBp] = c;
				}
			}
		}
		return refSeq;
	}

	@Override
	public String getName() {
		return "Reads";
	}
}
