package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Exon;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Gene;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Transcript;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.Strand;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;
import fi.csc.microarray.constants.VisualConstants;

/**
 * Track for showing transcripts. Lower zoom level version of {@link GeneTrack}.
 *
 */
public class TranscriptTrack extends Track {

	private Collection<Gene> genes = new TreeSet<Gene>();

	List<Integer> occupiedSpace = new ArrayList<Integer>();

	private Color color;

	public enum PartColor {
		CDS(VisualConstants.COLOR_BLUE), UTR(VisualConstants.COLOR_ORANGE), START_CODON(Color.gray);

		//		CDS(new Color(64, 192, 64)), UTR(new Color(192, 64, 64)), START_CODON(Color.gray);
		public Color c;

		PartColor(Color c) {
			this.c = c;
		}
	}

	public TranscriptTrack(View view, DataSource file, Color color, long maxBpLength) {

		super(view, file);
		this.color = color;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		occupiedSpace.clear();

		if (genes != null) {
			
			List<Drawable> geneDrawables = new ArrayList<Drawable>();

			//Use iterator to be able to remove genes that are out of sight
			Iterator<Gene> geneIter = genes.iterator();			
			Gene gene;

			while (geneIter.hasNext()) {
				gene = geneIter.next();
				
				if (!gene.getRegion().intersects(getView().getBpRegion())) {
					
					geneIter.remove();
					continue;
				}

				//Transcript collection refers to original data from the data layer, so out-of-sight
				//transcripts can't be removed
				for (Transcript transcript : gene.getTranscripts()) {		
					
					Rectangle rect = new Rectangle();

					rect.x = getView().bpToTrack(transcript.getRegion().start);
					int x = rect.x;
					rect.width = getView().bpToTrack(transcript.getRegion().end) - rect.x;
					int x2 = getView().bpToTrack(transcript.getRegion().end);

					int i = 0;

					while (occupiedSpace.size() > i && occupiedSpace.get(i) > rect.x) {
						i++;
					}

					int end = rect.x + rect.width;

					if (occupiedSpace.size() > i) {
						occupiedSpace.set(i, end);
					} else {
						occupiedSpace.add(end);
					}

					rect.y = (int) (((i + 1) * (14)));
					int y = rect.y + 2;
					rect.height = 2;

					drawables.add(new LineDrawable(x, y, x2, y, Color.darkGray));

					rect.height = 4;

					// draw arrow
					if (transcript.getRegion().getStrand() == Strand.REVERSED) {
						drawables.addAll(getArrowDrawables(rect.x, rect.y, -rect.height, rect.height));
					} else {
						drawables.addAll(getArrowDrawables(rect.x + rect.width, rect.y, rect.height, rect.height));
					}

					String name = transcript.getName();

					if (isNameVisible(rect)) {
						
						if (name == null) {
							name = "n/a";
						}
						
						drawTextAboveRectangle(name, drawables, rect, 1);
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
							c = PartColor.START_CODON.c;
							break;
						case STOP_CODON:
							// TODO Check how this should be visualised
							c = PartColor.UTR.c;
							break;
						default:
							System.err.println("Gene description not recognised: " + feature);
							c = Color.blue;
						}

						rect.x = getView().bpToTrack(exon.getRegion().start);
						rect.width = getView().bpToTrack(exon.getRegion().end) - rect.x;
						rect.height = 4;

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

	public void processAreaResult(AreaResult areaResult) {

		for (RegionContent content : areaResult.getContents()) {
			
			// Sorting is needed to draw partly overlapping genes in the same order every time
			if (!areaResult.getStatus().concise && content.region.getStrand() == getStrand()) {

				Gene gene = (Gene) content.values.get(ColumnType.VALUE);
				
				genes.add(gene);

			}
		}
		getView().redraw();
	}

	private long maxBpLength;

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
		// hide if visible region is too large
		return (super.isVisible() &&
				getView().getBpRegion().getLength() <= maxBpLength);
	}

	@Override
	public Map<DataSource, Set<ColumnType>> requestedData() {
		HashMap<DataSource, Set<ColumnType>> datas = new
				HashMap<DataSource, Set<ColumnType>>();
		datas.put(file, new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {
				ColumnType.VALUE })));
		return datas;
	}

	@Override
	public boolean isConcised() {
		return false;
	}
}
