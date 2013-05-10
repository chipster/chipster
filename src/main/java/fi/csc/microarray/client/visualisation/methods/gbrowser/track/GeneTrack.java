package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.TreeMap;

import fi.csc.microarray.client.visualisation.methods.gbrowser.fileIndex.GtfToFeatureConversion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.LayoutTool.LayoutMode;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Exon;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Gene;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.GeneSet;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.PositionAndStringKey;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Track for genes. Higher zoom level version of {@link TranscriptTrack}.
 *
 */
public class GeneTrack extends Track {

	private HashSet<Exon> exons = new HashSet<Exon>();
	private List<Integer> occupiedSpace = new ArrayList<Integer>();

	private Color color;


	public GeneTrack(Color color) {

		this.color = color;
		this.layoutMode = this.defaultLayoutMode = LayoutMode.FULL;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		occupiedSpace.clear();
				
		if (exons != null) {
			
			GeneSet geneSet = new GeneSet();				
			geneSet.add(exons.iterator(), view.getRequestRegion().grow(GtfToFeatureConversion.MAX_INTRON_LENGTH * 2));		
			
			TreeMap<PositionAndStringKey, Gene> sortedGenes = new TreeMap<PositionAndStringKey, Gene>();

			Iterator<Gene> geneIter = geneSet.values().iterator();
			while (geneIter.hasNext()) {

				Gene gene = geneIter.next();

				if (!getView().requestIntersects(gene.getRegion())) {
					geneIter.remove();
					continue;
				}

				PositionAndStringKey key = new PositionAndStringKey(gene.getRegion().start, gene.getId());
				sortedGenes.put(key, gene);
			}

			for (Gene gene : sortedGenes.values()) {
				
				if (!getView().getBpRegion().intersects(gene.getRegion())) {
					continue;
				}
				
				String name = null;
				
				if (gene.getName() != null) {
					name = gene.getName();
				} else if (gene.getId() != null) {
					name = gene.getId();
				} else {
					name = "n/a";
				}

				createDrawable(gene.getRegion().start, gene.getRegion().end, 10, color, name, drawables);
			}
		}
		
		return drawables;
	}

	private void createDrawable(BpCoord startBp, BpCoord endBp, int height, Color c, String name, Collection<Drawable> drawables) {
		Rectangle rect = new Rectangle();

		rect.x = getView().bpToTrack(startBp);
		rect.width = getView().bpToTrack(endBp) - rect.x;

		int i = 0;

		while (occupiedSpace.size() > i && occupiedSpace.get(i) > rect.x + 1) {
			i++;
		}

		int end = rect.x + rect.width;

		if (occupiedSpace.size() > i) {
			occupiedSpace.set(i, end);
		} else {
			occupiedSpace.add(end);
		}

		rect.y = (int) ((i + 1) * (height + 2));
		rect.height = height;

		drawables.add(new RectDrawable(rect, c, null));
		if (isNameVisible(rect)) {

			// draw name to leftmost visible part of the gene rectangle
			drawTextAboveRectangle(name, drawables, rect, 10);
		}
	}

	public void processDataResult(DataResult dataResult) {

		for (RegionContent content : dataResult.getContents()) {

			Object value = content.values.get(DataType.VALUE);
			
			if (value instanceof Exon) {
				Exon exon = (Exon)value;

				if (exon.getRegion().getStrand() == getStrand()) {

					this.exons.add(exon);
				}
			}
		}
	}
    
    @Override
	public void defineDataTypes() {
		addDataType(DataType.VALUE);
	}
    	
	@Override
	public int getMinHeight() {
		return 100;
	}
}
