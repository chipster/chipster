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

import fi.csc.microarray.client.visualisation.methods.gbrowser.ChunkDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Track for genes.
 *
 */
public class GeneTrack extends Track {

	private Collection<RegionContent> reads = new TreeSet<RegionContent>();
	private List<Integer> occupiedSpace = new ArrayList<Integer>();

	private long maxBpLength;
	private long minBpLength;

	private boolean wasLastConcised = true;
	private Color color;


	public GeneTrack(View view, ChunkDataSource file, Class<? extends AreaRequestHandler> handler,
	        Color color, long minBpLength, long maxBpLength) {
		super(view, file, handler);
		this.color = color;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		occupiedSpace.clear();

		if (reads != null) {

			Iterator<RegionContent> iter = reads.iterator();
			while (iter.hasNext()) {

				RegionContent read = iter.next();

				// FIXME this and all the other incarnations of the same 3 lines should be refactored up to Track or something
				if (!read.region.intercepts(getView().getBpRegion())) {
					iter.remove();
					continue;
				}

				String name = ((String) read.values.get(ColumnType.DESCRIPTION));

				createDrawable(read.region.start, read.region.end, 10, color, name, drawables);
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

		rect.y = (int) (getView().getTrackHeight() - ((i + 1) * (height + 2)));
		rect.height = height;

		drawables.add(new RectDrawable(rect, c, null));
		if (rect.width > name.length() * 7) {

			// TODO fix the extra quote mark in file
			name = name.replaceAll("\"", "");

			drawables.add(new TextDrawable(rect.x, rect.y + 10, name, Color.DARK_GRAY));
		}
	}

	public void processAreaResult(AreaResult<RegionContent> areaResult) {

		if (areaResult.status.concise == this.isConcised() && areaResult.content.values.get(ColumnType.STRAND) == getStrand()) {

			this.reads.add(areaResult.content);

			getView().redraw();
		}
	}


	@Override
	public void updateData() {
		if (wasLastConcised != isConcised()) {
			reads.clear();
			wasLastConcised = isConcised();
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
    public Map<DataSource, Set<ColumnType>> requestedData() {
        HashMap<DataSource, Set<ColumnType>> datas = new
        HashMap<DataSource, Set<ColumnType>>();
        datas.put(file, new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {
                ColumnType.STRAND,
                ColumnType.DESCRIPTION,
                ColumnType.VALUE })));
        return datas;
    }

	@Override
	public boolean isConcised() {
		return false;
	}
}
