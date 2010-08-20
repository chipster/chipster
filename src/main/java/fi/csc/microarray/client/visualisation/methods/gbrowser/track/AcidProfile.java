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
public class AcidProfile extends Track {

	private long maxBpLength;
	private long minBpLength;
	
	private enum Acid {
		A, C, G, T 
	}
	
	private Collection<RegionContent> forwardReads = new TreeSet<RegionContent>();
	
	static private Map<Acid, Color> colorMap;
	
	static {
		
		colorMap = new HashMap<Acid, Color>();

		colorMap.put(Acid.A, new Color(64, 192, 64));
		colorMap.put(Acid.C, new Color(64, 64, 192));
		colorMap.put(Acid.G, new Color(128, 128, 128));
		colorMap.put(Acid.T, new Color(192, 64, 64));
	}

	public AcidProfile(View view, DataSource file, Class<? extends AreaRequestHandler> handler,
			long minBpLength, long maxBpLength) {
		
		super(view, file, handler);

		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}
	
	/**
	 * Get drawables for a collection of reads.
	 * 
	 * @return
	 */
	private Collection<Drawable> getDrawableReads(Collection<RegionContent> reads) {
        
        TreeMap<Long, Long> ACollector = new TreeMap<Long, Long>();
        TreeMap<Long, Long> CCollector = new TreeMap<Long, Long>();
        TreeMap<Long, Long> GCollector = new TreeMap<Long, Long>();
        TreeMap<Long, Long> TCollector = new TreeMap<Long, Long>();
        
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
            
            
            
			String seq = ((String) read.values.get(ColumnType.SEQUENCE));
			
			int seqLength = seq.length();

			if (seq != null) {
				seq = seq.trim().toUpperCase();			

				for (Long i = 0l; i < seqLength; i++) {
					
					Long j = i + read.region.start.bp;
					
					TreeMap<Long, Long> collector = null;
					
					switch(seq.charAt(i.intValue())) {
						case 'A':
							collector = ACollector;
							break;
						case 'C':
							collector = CCollector;
							break;
						case 'G':
							collector = GCollector;
							break;
						case 'T':
							collector = TCollector;
							break;
					}

					if (collector != null) {
						if (collector.containsKey(j)) {
							collector.put(j, collector.get(j) + 1);
						} else {
							collector.put(j, 1L);
						}
					}
				}
			}
        }
        
        
        Collection<Drawable> drawables = getEmptyDrawCollection();
        
        createDrawables(ACollector, colorMap.get(Acid.A), drawables, lastChromosome);
        createDrawables(CCollector, colorMap.get(Acid.C), drawables, lastChromosome);
        createDrawables(GCollector, colorMap.get(Acid.G), drawables, lastChromosome);
        createDrawables(TCollector, colorMap.get(Acid.T), drawables, lastChromosome);
        
        return drawables;
        
	}
		
	private Collection<Drawable> createDrawables(TreeMap<Long, Long> collector, Color color, 
			Collection<Drawable> drawables, Chromosome lastChromosome) {
		
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
		drawables.addAll(getDrawableReads(forwardReads));

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
	               // implement later if necessary, bam converts anything to forward strand anyway
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
