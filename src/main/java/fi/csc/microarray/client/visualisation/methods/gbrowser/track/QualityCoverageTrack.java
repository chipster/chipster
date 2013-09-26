package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeMap;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Cigar;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;

/**
 * Track for showing the coverage of reads. Profile is drawn by calculating
 * the number of nucleotides hitting each basepair location. Should look
 * similar to IntensityTrack, but is exact.
 * 
 * If track's strand is set to {@link Strand#BOTH}, two profiles are drawn, one for
 * each strand.
 *
 * @see CoverageEstimateTrack
 */
public class QualityCoverageTrack extends Track {

	//Contains also reverse complemented reverse reads
	private Collection<Feature> forwardReads = new TreeSet<Feature>();
	private Color forwardColor;

	public QualityCoverageTrack(Color forwardColor) {
		super();

		this.forwardColor = forwardColor;
		
		this.setStrand(Strand.BOTH);
	}

	/**
	 * Get drawables for a collection of reads.
	 * 
	 * @return
	 */
	private Collection<Drawable> getDrawableReads(Collection<Feature> reads, Color color) {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		Chromosome chr = getView().getBpRegion().start.chr;

		TreeMap<Long, Float> collector = getQualities(reads);

		// width of a single bp in pixels
		int bpWidth = (int) (getView().getWidth() / getView().getBpRegion().getLength());

		// maximum y coordinate
		int bottomlineY = 0;		

		// prepare lines that make up the profile for drawing
		Iterator<Long> bpLocations = collector.keySet().iterator();
		if (bpLocations.hasNext()) {
			Long lastBpLocation = bpLocations.next();

			// draw a line from the beginning of the graph to the first location
			int startX = getView().bpToTrack(new BpCoord(lastBpLocation, chr));
			long startY = (long)(float)collector.get(lastBpLocation);
			drawables.add(new LineDrawable(0, bottomlineY,
					(int)(startX - bpWidth), bottomlineY, color));
			drawables.add(new LineDrawable((int)(startX - bpWidth), bottomlineY,
					startX, (int)(bottomlineY + startY), color));

			// draw lines for each bp region that has some items
			while (bpLocations.hasNext()) {
				Long currentBpLocation = bpLocations.next();

				startX = getView().bpToTrack(new BpCoord(lastBpLocation, chr));
				startY = (long)(float)collector.get(lastBpLocation);
				int endX = getView().bpToTrack(new BpCoord(currentBpLocation, chr));
				long endY = (long)(float)collector.get(currentBpLocation);

				// TODO could be approximated using natural cubic spline interpolation,
				//      then having a formula S(x) for each interval we could draw
				//      several lines approximating the S(x)

				if (currentBpLocation - lastBpLocation == 1) {
					// join adjacent bp locations with a line
					drawables.add(new LineDrawable(startX, (int)(bottomlineY + startY),
							endX, (int)(bottomlineY + endY), color));
				} else {
					// join locations that are more than one bp apart
					drawables.add(new LineDrawable((int)startX, (int)(bottomlineY + startY),
							(int)(startX + bpWidth), bottomlineY, color));
					drawables.add(new LineDrawable((int)(startX + bpWidth), bottomlineY,
							(int)(endX - bpWidth), bottomlineY, color));
					drawables.add(new LineDrawable((int)(endX - bpWidth), bottomlineY,
							(int)endX, (int)(bottomlineY + endY), color));
				}

				lastBpLocation = currentBpLocation;
			}

			// draw a line from the last location to the end of the graph
			int endX = getView().bpToTrack(new BpCoord(lastBpLocation, chr));
			long endY = (long)(float)collector.get(lastBpLocation);
			drawables.add(new LineDrawable(endX, (int)(bottomlineY + endY),
					(int)(endX + bpWidth), bottomlineY, color));
			drawables.add(new LineDrawable((int)(endX + bpWidth), bottomlineY,
					getView().getWidth(), bottomlineY, color));
		}

		collector.clear();

		return drawables;
	}
	
	private TreeMap<Long, Float> getQualities(Collection<Feature> reads) {

		TreeMap<Long, Float> collector = new TreeMap<Long, Float>();
		Iterator<Feature> iter = reads.iterator();

		// iterate over RegionContent objects (one object corresponds to one read)
		while (iter.hasNext()) {

			Feature read = iter.next();

			// remove those that are not in this view
			if (!read.region.intersects(getView().getBpRegion())) {
				iter.remove();
				continue;
			}

			Cigar cigar = (Cigar) read.values.get(DataType.CIGAR);

			String quality = ((String) read.values.get(DataType.QUALITY));
			
			if (quality == null) {
				continue;
			}

			for (int i = 0; i < quality.length(); i++) {

				long refIndex = cigar.getReferenceIndex(i);

				if (refIndex == -1) {
					//Skip insertions
					continue;
				}

				Long bp = refIndex + read.region.start.bp;
				
				char qualityChar = quality.charAt(i);
				
				//Scale to approximately same size than read coverage
				float qualityValue = ((int)qualityChar - 33) / 25f;

				if (!collector.containsKey(bp)) {
					collector.put(bp, qualityValue);
				} else {
					collector.put(bp, collector.get(bp) + qualityValue);
				}
			}
		}	        
		return collector;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		// add drawables from both reads (if present)
		drawables.addAll(getDrawableReads(forwardReads, forwardColor));

		return drawables;
	}

	public void processDataResult(DataResult dataResult) {

		for (Feature content : dataResult.getFeatures()) {

			// check that dataResult has 
			// correct strand
			if (getStrand() == content.values.get(DataType.STRAND) || 
					getStrand() == Strand.BOTH) {

				forwardReads.add(content);
			}
		}
	}

	@Override
	public int getTrackHeight() {
		return 100;
	}
	
    @Override
	public void defineDataTypes() {
		addDataType(DataType.ID);
		addDataType(DataType.STRAND);
		addDataType(DataType.QUALITY);
		addDataType(DataType.CIGAR);
	}

	/**
	 * @see GBrowserView#drawView
	 */
	@Override
	public boolean canExpandDrawables() {
		return true;
	}

	@Override
	public String getTrackName() {
		return "QualityCoverageTrack";
	}
}
