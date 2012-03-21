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

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.ReadPart;

/**
 * Track for showing the coverage of reads. Profile is drawn by calculating
 * the number of nucleotides hitting each basepair location. Should look
 * similar to IntensityTrack, but is exact.
 * 
 * If track's strand is set to {@link Strand#BOTH}, two profiles are drawn, one for
 * each strand.
 *
 * @see IntensityTrack
 */
public class CoverageTrack extends Track {

	private long maxBpLength;
	private long minBpLength;

	private Color forwardColor;
	private Color backwardColor;
	private ReadpartDataProvider readpartProvider;


	public CoverageTrack(View view, DataSource file, ReadpartDataProvider readpartProvider,
	        Color forwardColor, Color backwardColor, long minBpLength, long maxBpLength) {
		super(view, file);
		this.forwardColor = forwardColor;
		this.backwardColor = backwardColor;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
		this.readpartProvider = readpartProvider;
	}
	
	/**
	 * Get drawables for a collection of reads.
	 * 
	 * @return
	 */
	private Collection<Drawable> getDrawableReads(Strand strand, Color color) {
        Collection<Drawable> drawables = getEmptyDrawCollection();
        
        TreeMap<Long, Long> collector = new TreeMap<Long, Long>();
        Chromosome lastChromosome = null;

		Iterable<ReadPart> readParts = readpartProvider.getReadparts(strand); 

		for (ReadPart element : readParts) {

			// Skip elements that are not in this view
			if (!element.intersects(getView().getBpRegion())) {
				continue;
			}

			// collect relevant data for this read
			BpCoord startBp = element.start;
			BpCoord endBp = element.end;
			lastChromosome = element.start.chr;

			int seqLength = (int) (endBp.minus(startBp) + 1);

			for (Long i = element.start.bp; i <= (element.start.bp + seqLength); i++) {
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
        int baselineY = 0;

        // prepare lines that make up the profile for drawing
        Iterator<Long> bpLocations = collector.keySet().iterator();
        if (bpLocations.hasNext()) {
            Long lastBpLocation = bpLocations.next();
            
            // draw a line from the beginning of the graph to the first location
            int startX = getView().bpToTrack(new BpCoord(lastBpLocation, lastChromosome));
            long startY = collector.get(lastBpLocation);
            drawables.add(new LineDrawable(0, baselineY,
                    (int)(startX - bpWidth), baselineY, color));
            drawables.add(new LineDrawable((int)(startX - bpWidth), baselineY,
                    startX, (int)(baselineY + startY), color));

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
                    drawables.add(new LineDrawable(startX, (int)(baselineY + startY),
                            endX, (int)(baselineY + endY), color));
                } else {
                    // join locations that are more than one bp apart
                    drawables.add(new LineDrawable((int)startX, (int)(baselineY + startY),
                            (int)(startX + bpWidth), baselineY, color));
                    drawables.add(new LineDrawable((int)(startX + bpWidth), baselineY,
                            (int)(endX - bpWidth), baselineY, color));
                    drawables.add(new LineDrawable((int)(endX - bpWidth), baselineY,
                            (int)endX, (int)(baselineY + endY), color));
                }
                
                lastBpLocation = currentBpLocation;
            }

            // draw a line from the last location to the end of the graph
            int endX = getView().bpToTrack(new BpCoord(lastBpLocation, lastChromosome));
            long endY = collector.get(lastBpLocation);
            drawables.add(new LineDrawable(endX, (int)(baselineY + endY),
                    (int)(endX + bpWidth), baselineY, color));
            drawables.add(new LineDrawable((int)(endX + bpWidth), baselineY,
                    getView().getWidth(), baselineY, color));
        }
        
        collector.clear();
        
        return drawables;
	}

	@Override
	public Collection<Drawable> getDrawables() {
        Collection<Drawable> drawables = getEmptyDrawCollection();

        // add drawables from both reads (if present)
		drawables.addAll(getDrawableReads(Strand.FORWARD, forwardColor));
		drawables.addAll(getDrawableReads(Strand.REVERSED, backwardColor));

		return drawables;
	}

	public void processAreaResult(AreaResult areaResult) {
		// Do not listen to actual read data, because that is taken care by ReadpartDataProvider
	}

    @Override
    public Integer getHeight() {
        if (isVisible()) {
            // return super.getHeight();
        	return 100;
        } else {
            return 0;
        }
    }
    
    @Override
    public boolean isStretchable() {
    	
    	return false;
        // stretchable unless hidden
        // return isVisible();
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
	
	@Override
	public String getName() {
		return "ProfileTrack";
	}
}
