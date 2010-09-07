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
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.client.visualisation.methods.gbrowser.utils.Sequence;

/**
 * Track that shows actual content of reads. At high level, values of nucleotides are not shown. At low level, they are
 * shown.
 *
 */
public class SeqBlockTrack extends Track {

	public static final String DUMMY_SEQUENCE = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
								+ "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
								+ "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";

	private static final int SPACE_BETWEEN_READS = 2;

	private static final Color[] charColors = new Color[] {
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

	public SeqBlockTrack(View view, DataSource file, 
			Class<? extends AreaRequestHandler> handler,
			Color fontColor, long minBpLength,
			long maxBpLength) {
		super(view, file, handler);
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		occupiedSpace.clear();
		
		// if SNP highlight mode is on, we need reference sequence data
		char[] refSeq = getReferenceArray();

		// iterate over RegionContent objects (one object corresponds to one read)
		Iterator<RegionContent> iter = reads.iterator();

		while (iter.hasNext()) {

			RegionContent read = iter.next();

			// remove those that are not in this view
			if (!read.region.intercepts(getView().getBpRegion())) {
				iter.remove();
				continue;
			}

			// collect relevant data for this read
			Color blockColor = Color.gray;
			int height = 4;
			BpCoord startBp = read.region.start;
			BpCoord endBp = read.region.end;
			String seq = ((String) read.values.get(ColumnType.SEQUENCE));

			if (seq != null) {
				// we have sequence
				seq = seq.trim();
				
			} else {
				// we don't have a sequence, generate something
				int seqLength = (int) (endBp.minus(startBp) + 1);
				
				if (seqLength > 128) {
					// TODO what happened?
					seqLength = 0;
				}
				
				seq = DUMMY_SEQUENCE.substring(0, seqLength);
			}
			
            // width in bp, sequence length not always equal to endBp - startBp + 1
            int widthBp = seq.length();

			// create rectangle covering the correct screen area (x-axis)
			Rectangle rect = new Rectangle();
			rect.x = getView().bpToTrack(startBp);
			rect.width = (int) Math.round(getView().bpWidth() * widthBp);

			// reads are drawn in order and placed in layers
			int layer = 0;
			while (occupiedSpace.size() > layer
					&& occupiedSpace.get(layer) > rect.x + 1) {
				layer++;
			}

			// read reserves the space of the layer from end to left corner of the screen
			int end = rect.x + rect.width;
			if (occupiedSpace.size() > layer) {
				occupiedSpace.set(layer, end);
			} else {
				occupiedSpace.add(end);
			}

			
			// now we can decide the y-axis
			rect.y = getYCoord(layer, height);
			rect.height = height;
			
			boolean lastBeforeMaxStackingDepthCut = getYCoord(layer + 1, height) < 0;
			
			// check maximum stacking depth
			if (rect.y < 0) {
				continue;
			}

			// reverse the read if on reverse strand
			if ((Strand) read.values.get(ColumnType.STRAND) == Strand.REVERSED) {
				StringBuffer buf = new StringBuffer(seq).reverse();
				seq = buf.toString();
			}

			
			if (rect.width < seq.length()) {
				// too little space - only show one rectangle for each read
			    
				Color color = blockColor;
				// mark last line that will be drawn
				if (lastBeforeMaxStackingDepthCut) {
					color = color.darker();
				}

				// no space to draw actual sequence
				drawables.add(new RectDrawable(rect, color, null));

			} else { 
				// enough space to show the actual sequence

				// draw arrow
				if (read.values.get(ColumnType.STRAND) == Strand.REVERSED) {
					drawables.addAll(getArrowDrawables(rect.x, rect.y,
							-rect.height, rect.height));
				} else {
					drawables.addAll(getArrowDrawables(rect.x + rect.width,
							rect.y, rect.height, rect.height));
				}

				// prepare to draw single nucleotides
				float x = getView().bpToTrackFloat(startBp);
				float increment = getView().bpWidth();
				float nextX;

				// draw something for each nucleotide
				for (int j = 0; j < seq.length(); j++) {

					char letter = seq.charAt(j);

					// character drawing commented away:
					//
					// if (rect.width > seq.length() * CHAR_WIDTH) {
					//
					// drawables.add(new TextDrawable((int) x + 1, rect.y +
					// 10, "" + letter, fontColor));
					//
					// // drawables.add(new TextDrawable((int)x, rect.y - 8,
					// // "" + letter, Color.black));
					//
					// }
					
					// choose a color depending on viewing mode
					Color bg = Color.white;
					int posInRef = startBp.bp.intValue() + j - getView().getBpRegion().start.bp.intValue();
					if (highlightSNP &&
					    posInRef >= 0 &&
					    posInRef < refSeq.length &&
					    Character.toLowerCase(refSeq[posInRef]) == Character.toLowerCase(letter)) {
					    bg = Color.gray;
					} else {
    					if (letter == 'A') {
    						bg = charColors[0];
    					} else if (letter == 'C') {
    						bg = charColors[1];
    					} else if (letter == 'G') {
    						bg = charColors[2];
    					} else if (letter == 'T') {
    						bg = charColors[3];
    					}
					}
					
					// mark last line that will be drawn
					if (lastBeforeMaxStackingDepthCut) {
						bg = bg.brighter();
					}

					nextX = x + increment;
					drawables.add(new RectDrawable(Math.round(x), rect.y,
					        Math.round(nextX) - Math.round(x), height, bg, null));
					x = nextX;
				}
			}
		}
		

		return drawables;
	}
	
	private int getYCoord(int layer, int height) {
		return (int) (getView().getTrackHeight() - ((layer + 1) * (height + SPACE_BETWEEN_READS)));
	}

	public void processAreaResult(AreaResult<RegionContent> areaResult) {

		// check that areaResult has same concised status (currently always false)
		// and correct strand
		if (areaResult.status.file == file &&
		    areaResult.status.concise == isConcised() &&
			areaResult.content.values.get(ColumnType.STRAND) == getStrand()) {
			
			// add this to queue of RegionContents to be processed
			this.reads.add(areaResult.content);
			getView().redraw();
		}
		
		// reference sequence data
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
        // stretchable unless hidden
        return isVisible();
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
                ColumnType.STRAND })));
        
        // We might also need reference sequence data
        if (highlightSNP) {
            datas.put(refData, new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {
                    ColumnType.SEQUENCE })));
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
     * Get reference sequence as an array, if set.
     */
    private char[] getReferenceArray() {
        char[] refSeq = new char[0];
        if (highlightSNP) {
            Iterator<RegionContent> iter = refReads.iterator();
            refSeq = new char[getView().getBpRegion().getLength().intValue() + 1];
            int startBp = getView().getBpRegion().start.bp.intValue();
            int endBp = getView().getBpRegion().end.bp.intValue();
            RegionContent read;
            while (iter.hasNext()) {
                read = iter.next();
                if (!read.region.intercepts(getView().getBpRegion())) {
                    iter.remove();
                    continue;
                }
                
                // we might need to reverse reference sequence
                char[] readBases;
                if (strand == Strand.REVERSED) {
                    readBases = Sequence.complement((String)read.values.get(ColumnType.SEQUENCE)).toCharArray();
                } else {
                    readBases = ((String)read.values.get(ColumnType.SEQUENCE)).toCharArray();
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
}
