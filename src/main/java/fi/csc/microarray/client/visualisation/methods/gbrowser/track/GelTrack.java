package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
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
import fi.csc.microarray.client.visualisation.methods.gbrowser.GenomePlot.ReadScale;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Gel-like track that shows intensity of some parameter as a continuous stripe
 * of gray shades or a heatmap. Not to be confused with Intensity track which
 * shows more something like distribution graph.
 * 
 * Basically this track shows the same information as ProfileTrack, but in a
 * different format.
 * 
 * If track's strand is set to Strand.BOTH, number of reads on both strands
 * are summed up.
 * 
 * @author naktinis
 *
 */
public class GelTrack extends Track {

    private long maxBpLength;
    private long minBpLength;

    private Collection<RegionContent> reads = new TreeSet<RegionContent>();
    private TreeMap<Long, Long> collector = new TreeMap<Long, Long>();
    private Color color;
    
    private Color BACKGROUND = Color.WHITE;

    public GelTrack(View view, DataSource file, Class<? extends AreaRequestHandler> handler,
            Color color, long minBpLength, long maxBpLength) {
        super(view, file, handler);
        this.color = color;
        this.minBpLength = minBpLength;
        this.maxBpLength = maxBpLength;
    }

    @Override
    public Collection<Drawable> getDrawables() {
        Collection<Drawable> drawables = getEmptyDrawCollection();

        collector.clear();
        
        // draw a black rectangle as the background
        drawables.add(new RectDrawable(new Rectangle(0, 0,
                getView().getWidth(), this.getHeight()), BACKGROUND, BACKGROUND));

        // iterate over RegionContent objects (one object corresponds to one read)
        if (reads != null) {
            Iterator<RegionContent> iter = reads.iterator();
            Chromosome lastChromosome = null;
            Long maxItems = 1L;
            
            // count items in each position
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
                        maxItems = Math.max(maxItems, collector.get(i) + 1);
                        collector.put(i, collector.get(i) + 1);
                    } else {
                        collector.put(i, 1L);
                    }
                }
            }            

            // prepare lines that make up the profile for drawing
            Iterator<Long> bpLocations = collector.keySet().iterator();
            if (bpLocations.hasNext()) {
                Long lastBpLocation = bpLocations.next();
                
                // hue
                float hue = Color.RGBtoHSB(color.getRed(), color.getGreen(),
                        color.getBlue(), null)[0];
                
                while (bpLocations.hasNext()) {
                    // take coordinates from bp
                    Long currentBpLocation = bpLocations.next();
                    long startX = getView().bpToTrack(new BpCoord(lastBpLocation, lastChromosome));
                    long endX = getView().bpToTrack(new BpCoord(currentBpLocation, lastChromosome));
                    
                    // choose lightness
                    Color c;
                    float lightness;
                    if (view.parentPlot.getReadScale() == ReadScale.AUTO) {
                        lightness = collector.get(currentBpLocation)/(float)maxItems; 
                    } else {
                        lightness = Math.min(collector.get(currentBpLocation)/
                                    (float)view.parentPlot.getReadScale().numReads, 1);
                    }
                    if (currentBpLocation - lastBpLocation == 1) {
                        c = Color.getHSBColor(hue, 0, 1-lightness);
                        
                        // draw a rectangle for each region
                        drawables.add(new RectDrawable(new Rectangle((int)startX, 0,
                                (int)(endX-startX), this.getHeight()), c, c));
                    }
                    
                    lastBpLocation = currentBpLocation;
                }
            }
        }

        return drawables;
    }

    public void processAreaResult(AreaResult<RegionContent> areaResult) {

        // check that areaResult has same concised status (currently always false)
        // and correct strand
        if (areaResult.status.concise == isConcised()
                && (getStrand() == areaResult.content.values.get(ColumnType.STRAND) ||
                    getStrand() == Strand.BOTH)) {
            
            // add this to queue of RegionContents to be processed
            this.reads.add(areaResult.content);
            getView().redraw();
        }
    }

    @Override
    public Integer getHeight() {
        if (isVisible()) {
            return 16;
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
    
    @Override
    public String getName() {
    	return "GelTrack";
    }
}
