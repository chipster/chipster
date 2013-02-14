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

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataSource.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutTool.LayoutMode;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Cigar;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.CigarItem.CigarItemType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ReadPart;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.util.Sequence;

/**
 * The read track, most important of all tracks. Shows actual content of reads using color coding.
 * 
 */
public class ReadPileTrack extends Track {
	
	public static final Color CUTOFF_COLOR = Color.ORANGE;
	
	private long maxBpLength;
	private long minBpLength;

	private DataSource refData;
	private Collection<RegionContent> refReads = new TreeSet<RegionContent>();

	private boolean highlightSNP = false;

	private DataSource readData;
	private Collection<RegionContent> reads = new TreeSet<RegionContent>();


	public ReadPileTrack(GBrowserView view, DataSource file, ReadpartDataProvider readpartProvider, Color fontColor, 
			long minBpLength, long maxBpLength) {
		super(view, file);
		this.readData = file;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
		this.layoutMode = this.defaultLayoutMode = LayoutMode.FILL;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		// If SNP highlight mode is on, we need reference sequence data
		char[] refSeq = highlightSNP ? getReferenceArray(refReads, view, strand) : null;

		Iterator<RegionContent> readIter = reads.iterator();
		RegionContent read = null;
		
		TreeSet<RegionContent> dividedReads = new TreeSet<RegionContent>();
	
		Collection<CigarItemType> splitters = new HashSet<CigarItemType>();
		splitters.add(CigarItemType.N);

		while (readIter.hasNext()) {
			read = readIter.next();

			// Skip elements that are not in this view
			if (!read.region.intersects(getView().getBpRegion())) {
				readIter.remove();
				continue;
			}
			
			dividedReads.addAll(Cigar.splitRead(read, splitters));
		}

		// Preprocessing loop: Iterate over RegionContent objects (one object corresponds to one read)
		//Iterable<ReadPart> readParts = readpartProvider.getReadparts(getStrand()); 
		//int reads = 0;
		// Main loop: Iterate over ReadPart objects (one object corresponds to one continuous element)
		List<Integer> occupiedSpace = new ArrayList<Integer>();

		Iterator<RegionContent> splittedReadIter = dividedReads.iterator();
		RegionContent splittedRead = null;

		while (splittedReadIter.hasNext()) {
			splittedRead = splittedReadIter.next();

			// Skip elements that are not in this view
			if (!splittedRead.region.intersects(getView().getBpRegion())) {
				continue;
			}

			// Width in basepairs
			long widthInBps = splittedRead.region.getLength();

			// Create rectangle covering the correct screen area (x-axis)
			Rectangle readRect = new Rectangle();
			readRect.x = getView().bpToTrack(splittedRead.region.start);
			readRect.width = (int) Math.floor(getView().bpWidth() * widthInBps);

			// Do not draw invisible rectangles
			if (readRect.width < 2) {
				readRect.width = 2;
			}

			// Read parts are drawn in order and placed in layers
			int layer = 0;
			while (occupiedSpace.size() > layer && occupiedSpace.get(layer) > readRect.x + 1) {
				layer++;
			}

			// Read part reserves the space of the layer from end to left corner of the screen
			int end = readRect.x + readRect.width;
			if (occupiedSpace.size() > layer) {
				occupiedSpace.set(layer, end);
			} else {
				occupiedSpace.add(end);
			}
			
			// Now we can decide the y coordinate
			readRect.y = getYCoord(layer, GBrowserConstants.READ_HEIGHT);
			readRect.height = GBrowserConstants.READ_HEIGHT;

			boolean lastBeforeMaxStackingDepthCut = false;
			
			if (getLayoutMode() != LayoutMode.FULL) {
				// Check if we are about to go over the edge of the drawing area
				lastBeforeMaxStackingDepthCut = getYCoord(layer + 1, GBrowserConstants.READ_HEIGHT) > getHeight();
				
				// Check if we are over the edge of the drawing area
				if (readRect.y > getHeight()) {
					continue;
				}
			}
			

			for (ReadPart readPart : Cigar.splitElements(splittedRead)) {

				// Width in basepairs
				widthInBps = readPart.getLength();

				// Create rectangle covering the correct screen area (x-axis)
				Rectangle rect = new Rectangle();
				rect.x = getView().bpToTrack(readPart.start);
				rect.width = (int) Math.floor(getView().bpWidth() * widthInBps);
				rect.y = readRect.y;
				rect.height = readRect.height;

				// Do not draw invisible rectangles
				if (rect.width < 2) {
					rect.width = 2;
				}

				// Check if we have enough space for the actual sequence (at least pixel per nucleotide)
				String seq = readPart.getSequencePart();
				Cigar cigar = (Cigar) readPart.getRead().values.get(ColumnType.CIGAR);
				if (readPart.isVisible()) {
					if (rect.width < seq.length()) {
						// Too little space - only show one rectangle for each read part

						Color color = Color.gray;

						// Mark last line that will be drawn
						if (lastBeforeMaxStackingDepthCut) {
							color = CUTOFF_COLOR;
						}

						drawables.add(new RectDrawable(rect, color, color, cigar.toInfoString()));

					} else {
						// Enough space - show color coding for each nucleotide

						// Complement the read if on reverse strand
						if ((Strand) readPart.getRead().values.get(ColumnType.STRAND) == Strand.REVERSE) {

							StringBuffer buf = new StringBuffer(seq.toUpperCase());

							// Complement
							seq = buf.toString().replace('A', 'x'). // switch A and T
									replace('T', 'A').replace('x', 'T').

									replace('C', 'x'). // switch C and G
									replace('G', 'C').replace('x', 'G');
						}

						// Prepare to draw single nucleotides
						float increment = getView().bpWidth();
						float startX = getView().bpToTrackFloat(readPart.start);

						// Draw each nucleotide
						for (int j = 0; j < seq.length(); j++) {

							char letter = seq.charAt(j);

							long refIndex = j;

							// Choose a color depending on viewing mode
							Color bg = Color.white;
							Color border = Color.white;
							long posInRef = readPart.start.bp.intValue() + refIndex - getView().getBpRegion().start.bp.intValue();
							if (highlightSNP && posInRef >= 0 && posInRef < refSeq.length && Character.toLowerCase(refSeq[(int)posInRef]) == Character.toLowerCase(letter)) {
								bg = Color.gray;
								border = bg;
							} else {
								switch (letter) {
								case 'A':
									bg = GBrowserConstants.charColors[0];
									border = bg;
									break;
								case 'C':
									bg = GBrowserConstants.charColors[1];
									border = bg;
									break;
								case 'G':
									bg = GBrowserConstants.charColors[2];
									border = bg;
									break;
								case 'T':
									bg = GBrowserConstants.charColors[3];
									border = bg;
									break;
								case 'N':
									bg = Color.white;
									border = Color.gray;
									break;
								}
							}

							// Tell that we have reached max. stacking depth
							if (lastBeforeMaxStackingDepthCut) {
								bg = CUTOFF_COLOR;
								border = bg;
							}

							// Draw rectangle
							int x1 = Math.round(startX + ((float)refIndex) * increment);
							int x2 = Math.round(startX + ((float)refIndex + 1f) * increment);
							int width = Math.max(x2 - x1, 1);
							drawables.add(new RectDrawable(x1, rect.y, width, GBrowserConstants.READ_HEIGHT, bg, border, cigar.toInfoString()));

						}
					}
				} else if (readPart.getCigarItem() != null) { //invisible cigar type

					if (readPart.getCigarItem().getCigarItemType().equals(CigarItemType.D)) {
						Color color = Color.black;

						rect.grow(0, -1);

						drawables.add(new RectDrawable(rect, color, null, cigar.toInfoString()));
					}

					if (readPart.getCigarItem().getCigarItemType().equals(CigarItemType.I)) {
						Color color = Color.black;

						//drawables.add(new LineDrawable(rect.x, rect.y - 1, rect.x, rect.y + rect.height + 2, color));
						drawables.add(new RectDrawable(rect.x, rect.y - 2, 3, rect.height + 4, color, color));

					}					
				}
			}
		}

		return drawables;
	}

	private int getYCoord(int layer, int height) {
		return (int) ((layer + 1) * (height + GBrowserConstants.SPACE_BETWEEN_READS));
	}

	public void processAreaResult(AreaResult areaResult) {

		for (RegionContent regCont : areaResult.getContents()) {
			if (regCont.values.get(ColumnType.STRAND) == this.getStrand() && regCont.values.containsKey(ColumnType.SEQUENCE)) {
				this.reads.add(regCont);
			}
		}

		// "Spy" on reference sequence data, if available
		if (areaResult.getStatus().file == refData) {
			this.refReads.addAll(areaResult.getContents());
		}
	}

	@Override
	public boolean isVisible() {
		// visible region is not suitable
		return (super.isVisible() && getView().getBpRegion().getLength() > minBpLength && getView().getBpRegion().getLength() <= maxBpLength);
	}

	@Override
	public Map<DataSource, Set<ColumnType>> requestedData() {
		HashMap<DataSource, Set<ColumnType>> datas = new HashMap<DataSource, Set<ColumnType>>();
		datas.put(file, new HashSet<ColumnType>(Arrays.asList(new ColumnType[] { ColumnType.ID, ColumnType.SEQUENCE, ColumnType.STRAND, ColumnType.CIGAR })));

		// We might also need reference sequence data
		if (highlightSNP && this.getView().getBpRegion().getLength() < this.getView().getWidth() * 2) {
			datas.put(refData, new HashSet<ColumnType>(Arrays.asList(new ColumnType[] { ColumnType.SEQUENCE })));
		}

		return datas;
	}

	/**
	 * Enable SNP highlighting and set reference data.
	 * 
	 * @param highlightSNP
	 * @see ReadPileTrack.setReferenceSeq
	 */
	public void enableSNPHighlight(DataSource file) {
		// turn on highlighting mode
		highlightSNP = true;

		// set reference data
		refData = file;
		view.getQueueManager().addResultListener(file, this);
	}

	/**
	 * Disable SNP highlighting.
	 * 
	 * @param file
	 */
	public void disableSNPHiglight() {
		// turn off highlighting mode
		highlightSNP = false;
	}

	/**
	 * Convert reference sequence reads to a char array.
	 */
	public static char[] getReferenceArray(Collection<RegionContent> refReads, GBrowserView view, Strand strand) {
		char[] refSeq = new char[0];
		Iterator<RegionContent> iter = refReads.iterator();
		refSeq = new char[view.getBpRegion().getLength().intValue() + 1];
		int startBp = view.getBpRegion().start.bp.intValue();
		int endBp = view.getBpRegion().end.bp.intValue();
		RegionContent read;
		while (iter.hasNext()) {
			read = iter.next();
			if (!read.region.intersects(view.getBpRegion())) {
				iter.remove();
				continue;
			}

			// we might need to reverse reference sequence
			char[] readBases = null;
			if (read.values.get(ColumnType.SEQUENCE) != null) { //when showing negative coordinates
				if (strand == Strand.REVERSE) {
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
		}
		return refSeq;
	}

	@Override
	public String getName() {
		return "Reads";
	}
	
	@Override
	public int getMinHeight() {
		return 100;
	}
}
