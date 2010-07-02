package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Track for showing the coverage of reads. Profile is drawn by calculating
 * the number of nucleotides hitting each basepair location. Should look
 * similar to IntensityTrack, but is exact.
 *
 * @see IntensityTrack
 */
public class ProfileTrack extends Track {

	private long maxBpLength;
	private long minBpLength;

	private Collection<RegionContent> reads = new TreeSet<RegionContent>();
	private TreeMap<Long, Long> collector = new TreeMap<Long, Long>();
	private boolean wasLastConsised = true;
	private Color color;


	public ProfileTrack(View view, DataSource file, Class<? extends AreaRequestHandler> handler,
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

		// iterate over RegionContent objects (one object corresponds to one read)
		if (reads != null) {
			Iterator<RegionContent> iter = reads.iterator();
			Chromosome lastChromosome = null;
			
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

			// prepare lines that make up the profile for drawing
			Iterator<Long> bpLocations = collector.keySet().iterator();
			if (bpLocations.hasNext()) {
				Long lastBpLocation = bpLocations.next();
				
				// draw a line from the beginning of the graph to the first location
                int startX = getView().bpToTrack(new BpCoord(lastBpLocation, lastChromosome));
                long startY = collector.get(lastBpLocation);
                drawables.add(new LineDrawable(0, this.getHeight(),
                        (int)(startX - bpWidth), this.getHeight(), color));
                drawables.add(new LineDrawable((int)(startX - bpWidth), this.getHeight(),
                        startX, (int)(this.getHeight() - startY), color));

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
					    drawables.add(new LineDrawable(startX, (int)(this.getHeight() - startY),
					            endX, (int)(this.getHeight() - endY), color));
					} else {
					    // join locations that are more than one bp apart
                        drawables.add(new LineDrawable((int)startX, (int)(this.getHeight() - startY),
                                (int)(startX + bpWidth), this.getHeight(), color));
                        drawables.add(new LineDrawable((int)(startX + bpWidth), this.getHeight(),
                                (int)(endX - bpWidth), this.getHeight(), color));
                        drawables.add(new LineDrawable((int)(endX - bpWidth), this.getHeight(),
                                (int)endX, (int)(this.getHeight() - endY), color));
					}
					
					lastBpLocation = currentBpLocation;
				}

                // draw a line from the last location to the end of the graph
                int endX = getView().bpToTrack(new BpCoord(lastBpLocation, lastChromosome));
                long endY = collector.get(lastBpLocation);
                drawables.add(new LineDrawable(endX, (int)(this.getHeight() - endY),
                        (int)(endX + bpWidth), this.getHeight(), color));
                drawables.add(new LineDrawable((int)(endX + bpWidth), this.getHeight(),
                        getView().getWidth(), this.getHeight(), color));
			} else {
			    // there are no items, draw a straight line at zero
                drawables.add(new LineDrawable(0, this.getHeight(),
                        getView().getWidth(), this.getHeight(), color));
			}
		}

		return drawables;
	}

	public void processAreaResult(AreaResult<RegionContent> areaResult) {

		// check that areaResult has same concised status (currently always false)
		// and correct strand
		if (areaResult.status.concise == isConcised()
				&& areaResult.content.values.get(ColumnType.STRAND) == getStrand()) {
			
			// add this to queue of RegionContents to be processed
			this.reads.add(areaResult.content);
			getView().redraw();
		}
	}

	@Override
	public void updateData() {

		if (wasLastConsised != isConcised()) {
			reads.clear();
			wasLastConsised = isConcised();
		}
		super.updateData();
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
	public Collection<ColumnType> getDefaultContents() {
		return Arrays.asList(new ColumnType[] { ColumnType.SEQUENCE,
				ColumnType.STRAND });
	}

	@Override
	public boolean isConcised() {
		return false;
	}
	
	@Override
    public boolean canExpandDrawables() {
        return true;
    }
}
