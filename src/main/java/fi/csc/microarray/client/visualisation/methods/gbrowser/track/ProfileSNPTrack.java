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
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
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
 * If track's strand is set to Strand.BOTH, two profiles are drawn, one for
 * each strand.
 *
 * @see IntensityTrack
 */
public class ProfileSNPTrack extends Track {

	private long maxBpLength;
	private long minBpLength;

	private Collection<RegionContent> forwardReads = new TreeSet<RegionContent>();
	private Color forwardColor;

	private boolean highlightSNP = false;

	enum Acid { A, C, G, T };

	private Acid getAcid(char character) {
		switch(character) {
		case 'A':
			return Acid.A;
		case 'C':
			return Acid.C;
		case 'G':
			return Acid.G;
		case 'T':
			return Acid.T;
		}
		return null;
	}

	private class Base {

		public int[] acidCounts = new int[Acid.values().length];

		public Base() {
			Arrays.fill(acidCounts, 0);
		}

		public int getCoverage() {

			int sum = 0;

			for (Acid acid : Acid.values()) {
				sum += acidCounts[acid.ordinal()]; 
			}
			
			return sum;
		}
		
		public boolean isSNP() {
			//TODO Values should be compared to reference sequence, now we don't notice, if
			//all reads have the same mismatch
			
			int zeroCount = 0;
			
			for (Acid acid : Acid.values()) {
				if (acidCounts[acid.ordinal()] == 0) {
					zeroCount++;
				}
			}
			
			return zeroCount < 3;
		}
	}


	public ProfileSNPTrack(View view, DataSource file, Class<? extends AreaRequestHandler> handler,
			Color forwardColor, long minBpLength, long maxBpLength) {
		super(view, file, handler);
		this.forwardColor = forwardColor;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;

		setStrand(Strand.BOTH);
	}

	/**
	 * Get drawables for a collection of reads.
	 * 
	 * @return
	 */
	private Collection<Drawable> getDrawableReads(Collection<RegionContent> reads, Color color) {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		Chromosome chr = getView().getBpRegion().start.chr;

		// Count acids for each location
		TreeMap<Long, Base> collector = getAcidCounts(reads); 

		// Count width of a single bp in pixels
		int bpWidth = (int) (getView().getWidth() / getView().getBpRegion().getLength());

		// Count maximum y coordinate
		int maxY = this.getHeight() - 1;

		// prepare lines that make up the profile for drawing
		Iterator<Long> bpLocations = collector.keySet().iterator();
		if (bpLocations.hasNext()) {
			Long lastBpLocation = bpLocations.next();

			// draw a line from the beginning of the graph to the first location
			int startX = getView().bpToTrack(new BpCoord(lastBpLocation, chr));
			long startY = collector.get(lastBpLocation).getCoverage();
			drawables.add(new LineDrawable(0, maxY,
					(int)(startX - bpWidth), maxY, color));
			drawables.add(new LineDrawable((int)(startX - bpWidth), maxY,
					startX, (int)(maxY - startY), color));

			// draw lines for each bp region that has some items
			while (bpLocations.hasNext()) {
				Long currentBpLocation = bpLocations.next();
				
				startX = getView().bpToTrack(new BpCoord(lastBpLocation, chr));
				startY = collector.get(lastBpLocation).getCoverage();
				int endX = getView().bpToTrack(new BpCoord(currentBpLocation, chr));
				long endY = collector.get(currentBpLocation).getCoverage();

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
				
				Base base = collector.get(currentBpLocation);

				if (!highlightSNP || base.isSNP()) {
					int y = maxY;				

					for (Acid acid : Acid.values()) {

						int increment = base.acidCounts[acid.ordinal()];

						if (increment > 0) {
							Color c = SeqBlockTrack.charColors[acid.ordinal()];

							drawables.add(new RectDrawable(endX, y - increment, bpWidth, increment, c, c));

							y -= increment;
						}
					}
				}

				lastBpLocation = currentBpLocation;
			}

			// draw a line from the last location to the end of the graph
			int endX = getView().bpToTrack(new BpCoord(lastBpLocation, chr));
			long endY = collector.get(lastBpLocation).getCoverage();
			drawables.add(new LineDrawable(endX, (int)(maxY - endY),
					(int)(endX + bpWidth), maxY, color));
			drawables.add(new LineDrawable((int)(endX + bpWidth), maxY,
					getView().getWidth(), maxY, color));
		}

		collector.clear(); // FIXME don't clear, but keep'em cached (remember that there are two exit routes)

		return drawables;
	}

	/**
	 * Goes through data and gives count for each location and acid.
	 */
	private TreeMap<Long, Base> getAcidCounts(Collection<RegionContent> reads) {
	
		TreeMap<Long, Base> collector = new TreeMap<Long, Base>();
		Iterator<RegionContent> iter = reads.iterator();

		// iterate over RegionContent objects (one object corresponds to one read)
		while (iter.hasNext()) {

			RegionContent read = iter.next();

			// remove those that are not in this view
			if (!read.region.intercepts(getView().getBpRegion())) {
				iter.remove();
				continue;
			}

			String seq = (String) read.values.get(ColumnType.SEQUENCE);
			
			Cigar cigar = (Cigar) read.values.get(ColumnType.CIGAR);

			if (cigar != null) {
				for (int i = 0; i < seq.length(); i++) {

					Base base = null;

					long refIndex = cigar.getReferenceIndex(i);

					if (refIndex == -1) {
						//Skip insertions
						continue;
					}

					Long bp = refIndex + read.region.start.bp;

					if (!collector.containsKey(bp)) {
						base = new Base();
						collector.put(bp, base);
					} else {
						base = collector.get(bp);
					}

					Acid acid = getAcid(seq.charAt(i));

					if (acid != null) {
						base.acidCounts[acid.ordinal()]++;
					}
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

	public void processAreaResult(AreaResult<RegionContent> areaResult) {

		// check that areaResult has same concised status (currently always false)
		// and correct strand
		if (areaResult.status.file == file &&
				areaResult.status.concise == isConcised()) {

			// Don't care about strand
			forwardReads.add(areaResult.content);

			getView().redraw();
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
				ColumnType.STRAND,
				ColumnType.QUALITY,
				ColumnType.CIGAR})));
		
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
		return "ProfileSNPTrack";
	}

	public void enableSNPHighlight() {
		highlightSNP = true;
	}

	public void disableSNPHighlight() {
		highlightSNP = false;
	}
}
