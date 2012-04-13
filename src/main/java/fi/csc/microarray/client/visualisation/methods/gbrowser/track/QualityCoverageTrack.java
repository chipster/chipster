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
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Chromosome;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Cigar;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

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
public class QualityCoverageTrack extends Track {

	private long maxBpLength;
	private long minBpLength;

	//Contains also reverse complemented reverse reads
	private Collection<RegionContent> forwardReads = new TreeSet<RegionContent>();
	private Color forwardColor;

	public QualityCoverageTrack(View view, DataSource file, Color forwardColor, long minBpLength, long maxBpLength) {
		super(view, file);
		this.forwardColor = forwardColor;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
		
		this.setStrand(Strand.BOTH);
	}

	/**
	 * Get drawables for a collection of reads.
	 * 
	 * @return
	 */
	private Collection<Drawable> getDrawableReads(Collection<RegionContent> reads, Color color) {
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
	
	private TreeMap<Long, Float> getQualities(Collection<RegionContent> reads) {

		TreeMap<Long, Float> collector = new TreeMap<Long, Float>();
		Iterator<RegionContent> iter = reads.iterator();

		// iterate over RegionContent objects (one object corresponds to one read)
		while (iter.hasNext()) {

			RegionContent read = iter.next();

			// remove those that are not in this view
			if (!read.region.intersects(getView().getBpRegion())) {
				iter.remove();
				continue;
			}

			Cigar cigar = (Cigar) read.values.get(ColumnType.CIGAR);

			String quality = ((String) read.values.get(ColumnType.QUALITY));
			
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

	public void processAreaResult(AreaResult areaResult) {

		for (RegionContent content : areaResult.getContents()) {

			// check that areaResult has same concised status (currently always false)
			// and correct strand
			if (areaResult.getStatus().concise == isConcised()) {
				if (getStrand() == content.values.get(ColumnType.STRAND) || 
						getStrand() == Strand.BOTH) {

					forwardReads.add(content);
				}
			}
		}
		getView().redraw();
	}

	@Override
	public Integer getHeight() {
		if (isVisible()) {
			//return super.getHeight();
			
			return 100;
			
		} else {
			return 0;
		}
	}

	@Override
	public boolean isStretchable() {
		// stretchable unless hidden
		//return isVisible();
		
		return false;
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
				ColumnType.STRAND,
				ColumnType.QUALITY,
				ColumnType.CIGAR })));
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
		return "QualityCoverageTrack";
	}
}
