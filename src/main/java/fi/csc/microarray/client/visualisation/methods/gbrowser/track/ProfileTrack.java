package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeMap;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Track for showing the coverage of reads. Profile is drawn by calculating
 * the number of nucleotides hitting each basepair location. Should look
 * similar to IntensityTrack, but is exact.
 * 
 * If track's strand is set to Strand.BOTH, two profiles are drawn, one for
 * each strand.
 *
 * @see IntensityTrack
 */
public class ProfileTrack extends Track {

	private long maxBpLength;
	private long minBpLength;

	private Collection<RegionContent> forwardReads = new TreeSet<RegionContent>();
	private Collection<RegionContent> backwardReads = new TreeSet<RegionContent>();
	private Color forwardColor;
	private Color backwardColor;


	public ProfileTrack(View view, DataSource file, Class<? extends AreaRequestHandler> handler,
	        Color forwardColor, Color backwardColor, long minBpLength, long maxBpLength) {
		super(view, file, handler);
		this.forwardColor = forwardColor;
		this.backwardColor = backwardColor;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}
	
	/**
	 * Get drawables for a collection of reads.
	 * 
	 * @return
	 */
	private Collection<Drawable> getDrawableReads(Collection<RegionContent> reads, Color color) {
        Collection<Drawable> drawables = getEmptyDrawCollection();
        
        TreeMap<Long, Long> collector = new TreeMap<Long, Long>();
        Iterator<RegionContent> iter = reads.iterator();
        Chromosome lastChromosome = null;
        
        // iterate over RegionContent objects (one object corresponds to one read)
        while (iter.hasNext()) {

            RegionContent read = iter.next();

            // remove those that are not in this view
            if (!read.region.intercepts(getView().getBpRegion())) {
                iter.remove();
                continue;
            }

            // collect relevant data for this read
            BpCoord startBp = read.region.start;
            BpCoord endBp = read.region.end;
            lastChromosome = read.region.start.chr;
            
            int seqLength = (int) (endBp.minus(startBp) + 1);

            for (Long i = read.region.start.bp; i <= (read.region.start.bp + seqLength); i++) {
                if (collector.containsKey(i)) {
                    collector.put(i, collector.get(i) + 1);
                } else {
                    collector.put(i, 1L);
                }
            }
        }
        
        // width of a single bp in pixels
        int bpWidth = (int) (getView().getWidth() / getView().getBpRegion().getLength());
        
        // maximum y coordinate
        int maxY = this.getHeight() - 1;

        // prepare lines that make up the profile for drawing
        Iterator<Long> bpLocations = collector.keySet().iterator();
        if (bpLocations.hasNext()) {
            Long lastBpLocation = bpLocations.next();
            
            // draw a line from the beginning of the graph to the first location
            int startX = getView().bpToTrack(new BpCoord(lastBpLocation, lastChromosome));
            long startY = collector.get(lastBpLocation);
            drawables.add(new LineDrawable(0, maxY,
                    (int)(startX - bpWidth), maxY, color));
            drawables.add(new LineDrawable((int)(startX - bpWidth), maxY,
                    startX, (int)(maxY - startY), color));

            // draw lines for each bp region that has some items
            while (bpLocations.hasNext()) {
                Long currentBpLocation = bpLocations.next();

                startX = getView().bpToTrack(new BpCoord(lastBpLocation, lastChromosome));
                startY = collector.get(lastBpLocation);
                int endX = getView().bpToTrack(new BpCoord(currentBpLocation, lastChromosome));
                long endY = collector.get(currentBpLocation);
                   
                // TODO could be approximated using natural cubic spline interpolation,
                //      then having a formula S(x) for each interval we could draw
                //      several lines approximating the S(x)
                
                if (currentBpLocation - lastBpLocation == 1) {
                    // join adjacent bp locations with a line
                    drawables.add(new LineDrawable(startX, (int)(maxY - startY),
                            endX, (int)(maxY - endY), color));
                } else {
                    // join locations that are more than one bp apart
                    drawables.add(new LineDrawable((int)startX, (int)(maxY - startY),
                            (int)(startX + bpWidth), maxY, color));
                    drawables.add(new LineDrawable((int)(startX + bpWidth), maxY,
                            (int)(endX - bpWidth), maxY, color));
                    drawables.add(new LineDrawable((int)(endX - bpWidth), maxY,
                            (int)endX, (int)(maxY - endY), color));
                }
                
                lastBpLocation = currentBpLocation;
            }

            // draw a line from the last location to the end of the graph
            int endX = getView().bpToTrack(new BpCoord(lastBpLocation, lastChromosome));
            long endY = collector.get(lastBpLocation);
            drawables.add(new LineDrawable(endX, (int)(maxY - endY),
                    (int)(endX + bpWidth), maxY, color));
            drawables.add(new LineDrawable((int)(endX + bpWidth), maxY,
                    getView().getWidth(), maxY, color));
        }
        
        collector.clear();
        
        return drawables;
	}

	@Override
	public Collection<Drawable> getDrawables() {
        Collection<Drawable> drawables = getEmptyDrawCollection();

        // add drawables from both reads (if present)
		drawables.addAll(getDrawableReads(forwardReads, forwardColor));
		drawables.addAll(getDrawableReads(backwardReads, backwardColor));

		return drawables;
	}

	public void processAreaResult(AreaResult<RegionContent> areaResult) {

		// check that areaResult has same concised status (currently always false)
		// and correct strand
	    if (areaResult.status.concise == isConcised()) {
	        if (getStrand() == areaResult.content.values.get(ColumnType.STRAND) || 
	            getStrand() == Strand.BOTH) {
	            
	            // Put data to different collections for different strands
	            if (areaResult.content.values.get(ColumnType.STRAND) == Strand.FORWARD) {
	                forwardReads.add(areaResult.content);
	            } else {
	                backwardReads.add(areaResult.content);
	            }
	            getView().redraw();
	        }
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
}
