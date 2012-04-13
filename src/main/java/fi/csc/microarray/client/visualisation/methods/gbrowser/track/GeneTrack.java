package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.Gene;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Track for genes. Higher zoom level version of {@link TranscriptTrack}.
 *
 */
public class GeneTrack extends Track {

	private Collection<Gene> genes = new TreeSet<Gene>();
	private List<Integer> occupiedSpace = new ArrayList<Integer>();

	private long maxBpLength;
	private long minBpLength;

	private Color color;


	public GeneTrack(View view, DataSource file, Color color, long minBpLength, long maxBpLength) {
		super(view, file);
		this.color = color;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		occupiedSpace.clear();

		if (genes != null) {

			Iterator<Gene> iter = genes.iterator();
			while (iter.hasNext()) {

				Gene gene = iter.next();

				// FIXME this and all the other incarnations of the same 3 lines should be refactored up to Track or something
				if (!gene.getRegion().intersects(getView().getBpRegion())) {
					iter.remove();
					continue;
				}

				String name = gene.getName();

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
			
			if (name == null) {
				name = "n/a";
			}

			// draw name to leftmost visible part of the gene rectangle
			drawTextAboveRectangle(name, drawables, rect, 10);
		}
	}

	public void processAreaResult(AreaResult areaResult) {

		for (RegionContent content : areaResult.getContents()) {
			if (areaResult.getStatus().concise == this.isConcised()) {
				
				Gene gene = (Gene) content.values.get(ColumnType.VALUE);
								
				if (gene.getRegion().getStrand() == getStrand()) {
					this.genes.add(gene);

				}
			}			
		}
		getView().redraw();
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
                ColumnType.VALUE })));
        return datas;
    }

	@Override
	public boolean isConcised() {
		return false;
	}
}
