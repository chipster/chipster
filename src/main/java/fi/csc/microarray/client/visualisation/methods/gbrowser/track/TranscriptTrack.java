package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex.GtfToFeatureConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Exon;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Gene;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneSet;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.PositionAndStringKey;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Feature;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Transcript;

/**
 * Track for showing transcripts. Lower zoom level version of {@link GeneTrack}.
 *
 */
public class TranscriptTrack extends Track {

	private HashSet<Exon> exons = new HashSet<Exon>();

	List<Integer> occupiedSpace = new ArrayList<Integer>();

	public enum PartColor {
		CDS(GBrowserConstants.COLOR_BLUE), UTR(GBrowserConstants.COLOR_ORANGE);

		public Color c;

		PartColor(Color c) {
			this.c = c;
		}
	}

	public TranscriptTrack() {
		super();
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		occupiedSpace.clear();

		TreeMap<PositionAndStringKey, Transcript> sortedTranscripts = new TreeMap<PositionAndStringKey, Transcript>();

		if (exons != null) {
			
			GeneSet geneSet = new GeneSet();				
			geneSet.add(exons.iterator(), view.getRequestRegion().grow(GtfToFeatureConversion.MAX_INTRON_LENGTH * 2));						

			Iterator<Gene> iter = geneSet.values().iterator();
			while (iter.hasNext()) {
				
				//Use iterator to be able to remove genes that are out of sight
				Gene gene = iter.next();			

				if (!getView().requestIntersects(gene.getRegion())) {
					iter.remove();
					continue;
				}

				for (Transcript transcript : gene.getTranscripts()) {
					PositionAndStringKey key = new PositionAndStringKey(transcript.getRegion().start, transcript.getId());
					sortedTranscripts.put(key, transcript);
				}
			}

			List<Drawable> geneDrawables = new ArrayList<Drawable>();

			for (Transcript transcript : sortedTranscripts.values()) {
				
				if (!getView().getBpRegion().intersects(transcript.getRegion())) {
					continue;
				}

				Rectangle rect = new Rectangle();

				rect.x = getView().bpToTrack(transcript.getRegion().start);
				int x = rect.x;

				//End has to be increased by one, because the transcript includes the base at the end location
				BpCoord endCoord = new BpCoord(transcript.getRegion().end.bp + 1, transcript.getRegion().end.chr);
				rect.width = getView().bpToTrack(endCoord) - rect.x;
				int x2 = getView().bpToTrack(endCoord);

				int i = 0;

				while (occupiedSpace.size() > i && occupiedSpace.get(i) > rect.x) {
					i++;
				}

				int end = rect.x + rect.width;

				if (occupiedSpace.size() > i) {
					occupiedSpace.set(i, end + 1);
				} else {
					occupiedSpace.add(end + 1);
				}

				rect.y = (int) (((i + 1) * (14)));
				int y = rect.y + 2;
				rect.height = 2;

				drawables.add(new LineDrawable(x, y, x2, y, Color.darkGray));

				rect.height = 4;

				// draw arrow
				if (transcript.getRegion().getStrand() == Strand.REVERSE) {
					drawables.addAll(getArrowDrawables(rect.x, rect.y, -rect.height, rect.height));
				} else {
					drawables.addAll(getArrowDrawables(rect.x + rect.width, rect.y, rect.height, rect.height));
				}

				if (isNameVisible(rect)) {

					
					String name = null;
					
					if (transcript.getName() != null) {
						name = transcript.getName();
					} else if (transcript.getId() != null) {
						name = transcript.getId();
					} else {
						name = "n/a";
					}

					drawTextAboveRectangle(name, drawables, rect, 0);
				}

				for (Exon exon : transcript.getExons()) {

					//					if (part.values == null) {
					//						drawables.add(createDrawable(part.region.start, part.region.end, color));
					//					} else {					

					Exon.Feature feature = exon.getFeature();
					Color c;

					switch (feature) {
					case CDS:
						c = PartColor.CDS.c;
						break;
					case EXON:
						c = PartColor.UTR.c;
						break;
					case START_CODON:
						c = PartColor.CDS.c;
						break;
					case STOP_CODON:
						c = PartColor.CDS.c;
						break;
					case TRANSCRIPT:
						c = null;
						break;
					default:
						System.err.println("Gene description not recognised: " + feature);
						c = Color.gray;
					}

					
					rect.x = getView().bpToTrack(exon.getRegion().start);
					//End has to be increased by one, because the transcript includes the base at the end location
					BpCoord exonEnd = new BpCoord(exon.getRegion().end.bp + 1, exon.getRegion().end.chr);
					rect.width = getView().bpToTrack(exonEnd) - rect.x;
					rect.height = 4;

					if (c != null) {
						geneDrawables.add(new RectDrawable(rect, c, null));
					}
				}
			}

			Collections.sort(geneDrawables, new Comparator<Drawable>() {
				public int compare(Drawable one, Drawable other) {

					if (one.color.equals(PartColor.CDS.c) && other.color.equals(PartColor.UTR.c)) {
						return 1;
					} else if (one.color.equals(PartColor.UTR.c) && other.color.equals(PartColor.CDS.c)) {
						return -1;
					} else {
						return 0;
					}
				}
			});

			drawables.addAll(geneDrawables);

		}

		return drawables;
	}

	public void processDataResult(DataResult dataResult) {

		for (Feature content : dataResult.getFeatures()) {


				// Sorting is needed to draw partly overlapping genes in the same order every time
				if (content.region.getStrand() == getStrand()) {

					Object value = content.values.get(DataType.VALUE);
					
					if (value instanceof Exon) {
						Exon exon = (Exon)value;

					exons.add(exon);					
				}
			}
		}
	}
	
    @Override
	public void defineDataTypes() {
    	addDataType(DataType.VALUE);
	}
	
	@Override
	public int getTrackHeight() {
		return 100;
	}
	
	@Override
	public boolean isShowMoreCapable() {
		return true;
	}
}
