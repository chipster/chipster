package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutTool.LayoutMode;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Cigar;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.CigarItem.CigarItemType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ReadPart;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Sequence;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;

/**
 * The read track, most important of all tracks. Shows actual content of reads using color coding.
 * 
 */
public class ReadPileTrack extends Track {
	
	public static class Read extends Feature {
		
		private int layer;
		private boolean warning;

		public Read(Feature feature) {
			super(feature);
		}

		public boolean isWarning() {
			return warning;
		}

		public void setWarning(boolean warning) {
			this.warning = warning;
		}

		public int getLayer() {
			return layer;
		}

		public void setLayer(int layer) {
			this.layer = layer;
		}
	}
	
	public static final Color CUTOFF_COLOR = Color.ORANGE;

	private static final int MAX_FULL_HEIGHT = 1000;

	private DataThread refData;
	private Collection<Feature> referenceSequenceFeatures = new TreeSet<Feature>();

	private boolean highlightSNP = false;

	private TreeSet<Read> dividedReads = new TreeSet<>();

	private boolean markMultimappingReads;

	public ReadPileTrack(DataThread refData, Color fontColor) {
		super();
		this.refData = refData;
		this.layoutMode = this.defaultLayoutMode = LayoutMode.FILL;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		
		Collection<Drawable> drawables = getEmptyDrawCollection();

		// If SNP highlight mode is on, we need reference sequence data
		char[] refSeq = highlightSNP ? getReferenceArray(referenceSequenceFeatures, view, strand) : null;		
		

		Iterator<Read> splittedReadIter = dividedReads.iterator();
		Read splittedRead = null;
		
		while (splittedReadIter.hasNext()) {
			splittedRead = splittedReadIter.next();

			// Skip elements that are not visible
			if (!getView().viewIntersects(splittedRead.region)) {
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
			int layer = splittedRead.getLayer();
			boolean lastBeforeMaxStackingDepthCut = splittedRead.isWarning();
			
			// Now we can decide the y coordinate
			readRect.y = getYCoord(layer);
			readRect.height = GBrowserConstants.READ_HEIGHT;


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
				Cigar cigar = (Cigar) readPart.getRead().values.get(DataType.CIGAR);
				if (readPart.isVisible()) {
					if (rect.width < seq.length()) {
						// Too little space - only show one rectangle for each read part

						Color color = getDefaultReadColor(readPart.getRead());

						// Mark last line that will be drawn
						if (lastBeforeMaxStackingDepthCut) {
							color = CUTOFF_COLOR;
						}

						drawables.add(new RectDrawable(rect, color, color));

					} else {
						// Enough space - show color coding for each nucleotide

						// Complement the read if on reverse strand
						if ((Strand) readPart.getRead().values.get(DataType.STRAND) == Strand.REVERSE) {

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
								bg = getDefaultReadColor(readPart.getRead());
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

						drawables.add(new RectDrawable(rect, color, null));
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

	private Color getDefaultReadColor(Feature read) {
		if (markMultimappingReads) {
			if (read.values.containsKey(DataType.BAM_TAG_NH)) {
				Integer alignments = (Integer)read.values.get(DataType.BAM_TAG_NH);
				if (alignments != null && alignments > 1) {
					return Color.lightGray;
				}
			}
		}
		return Color.gray;
	}

	private int getYCoord(int layer) {
		return (int) ((layer + 1) * (GBrowserConstants.READ_HEIGHT + GBrowserConstants.SPACE_BETWEEN_READS));
	}
	
	private int getLayer(int yCoord) {
		return (int) (yCoord / (GBrowserConstants.READ_HEIGHT + GBrowserConstants.SPACE_BETWEEN_READS));
	}

	public void processDataResult(DataResult dataResult) {
		
		Collection<CigarItemType> splitters = new HashSet<CigarItemType>();
		splitters.add(CigarItemType.N);

		for (Feature feature : dataResult.getFeatures()) {
			if (feature.values.get(DataType.STRAND) == this.getStrand() && feature.values.containsKey(DataType.SEQUENCE)) {
				
				// Skip elements that are not in this view
				if (!getView().requestIntersects(feature.region)) {
					continue;
				}

				List<Feature> reads = Cigar.splitRead(feature, splitters);
				
				for (Feature featureRead : reads) {
					
					Read read = new Read(featureRead);
					dividedReads.add(read);
				}
				
			}
		}
		
		if (dataResult.getStatus().getDataThread() == refData) {			
			this.referenceSequenceFeatures.addAll(dataResult.getFeatures());
		}
		
		// do layout immediately to remove useless data
		doLayout();
	}
	
    private void doLayout() {
    	
		Iterator<Read> splittedReadIter = dividedReads.iterator();
		Read splittedRead = null;
		
		List<Long> occupiedSpace = new ArrayList<>();
		
		while (splittedReadIter.hasNext()) {
			splittedRead = splittedReadIter.next();

			if (!getView().requestIntersects(splittedRead.region)) {
				splittedReadIter.remove();
				continue;
			}

			// Read parts are drawn in order and placed in layers
			int layer = 0;
			while (occupiedSpace.size() > layer && occupiedSpace.get(layer) > splittedRead.region.start.bp + 1) {
				layer++;
			}
					
			// Read part reserves the space of the layer from end to left corner of the screen
			long end = splittedRead.region.end.bp;
			if (occupiedSpace.size() > layer) {
				occupiedSpace.set(layer, end);
			} else {
				occupiedSpace.add(end);
			}
									
			int maxHeight;
			
			if (getLayoutMode() == LayoutMode.FULL) {
								
				maxHeight = MAX_FULL_HEIGHT;				
			} else {
				maxHeight = getLayer(getLayoutComponent().getHeight());				
			}
			
			// Check if we are about to go over the edge of the drawing area
			boolean lastBeforeMaxStackingDepthCut = layer + 1 >= maxHeight;

			
			// Check if we are over the edge of the drawing area
			if (layer > maxHeight) {					
				splittedReadIter.remove();
				continue;
			}
			
			splittedRead.setLayer(layer);
			splittedRead.setWarning(lastBeforeMaxStackingDepthCut);	
		}
	}

	@Override
	public void defineDataTypes() {
		addDataType(DataType.ID);
		addDataType(DataType.SEQUENCE);
		addDataType(DataType.STRAND);
		addDataType(DataType.CIGAR);
		addDataType(DataType.BAM_TAG_NH);
		
		// We might also need reference sequence data
		if (highlightSNP && this.getView().getBpRegion().getLength() < this.getView().getWidth() * 2) {
			addDataType(DataType.SEQUENCE);
		}
	}

	/**
	 * Set SNP highlighting enabled or disabled.
	 * @param highlightSnp
	 * 
	 * @param highlightSNP
	 * @see ReadPileTrack.setReferenceSeq
	 */
	public void setSNPHighlight(boolean highlightSnp) {
		this.highlightSNP = highlightSnp;
	}

	/**
	 * Convert reference sequence reads to a char array.
	 */
	public static char[] getReferenceArray(Collection<Feature> refFeatures, GBrowserView view, Strand strand) {
		char[] refSeq = new char[0];
		Iterator<Feature> iter = refFeatures.iterator();
		refSeq = new char[view.getBpRegion().getLength().intValue() + 1];
		int startBp = view.getBpRegion().start.bp.intValue();
		int endBp = view.getBpRegion().end.bp.intValue();
		Feature feature;
		while (iter.hasNext()) {
			feature = iter.next();
			if (!view.requestIntersects(feature.region)) {
				iter.remove();
				continue;
			}

			// we might need to reverse reference sequence
			char[] readBases = null;
			if (feature.values.get(DataType.SEQUENCE) != null) { //when showing negative coordinates
				if (strand == Strand.REVERSE) {
					readBases = Sequence.complement((String) feature.values.get(DataType.SEQUENCE)).toCharArray();
				} else {
					readBases = ((String) feature.values.get(DataType.SEQUENCE)).toCharArray();
				}

				int readStart = feature.region.start.bp.intValue();
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
	public String getTrackName() {
		return "Reads";
	}
	
	@Override
	public int getTrackHeight() {
		return 100;
	}
	
	@Override	
	public void initializeListener() {
		super.initializeListener();
		
		// Add listener for reference file
		if (dataThreads != null && refData != null) {
			view.getQueueManager().addDataResultListener(refData, this);
		}
	}
	
	@Override
	public boolean isShowMoreCapable() {
		return true;
	}

	public void setMarkMultimappingReads(boolean markMultimappingReads) {
		this.markMultimappingReads = markMultimappingReads;
	}
}
