package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Rectangle;
import java.util.Arrays;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.DataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.TabixDataSource;
import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class RepeatMaskerTrack extends Track{

	private long minBpLength;
	private long maxBpLength;
	private Color color;

	private Collection<RegionContent> regions = new TreeSet<RegionContent>();

	public RepeatMaskerTrack(View view, TabixDataSource file, long minBpLength, long maxBpLength){

		super(view, file);
		this.color = Color.lightGray;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {

		Collection<Drawable> drawables = getEmptyDrawCollection();      

		if (regions != null) {

			Iterator<RegionContent> iter = regions.iterator();

			RegionContent regionContent = null;

			while (iter.hasNext()) {

				regionContent = iter.next();

				if (!regionContent.region.intersects(getView().getBpRegion())) {
					iter.remove();
					continue;
				}

				BpCoord startBp = regionContent.region.start;
				//Increase by one to fill the last base pair
				BpCoord endBp = new BpCoord(regionContent.region.end.bp + 1, regionContent.region.end.chr);
				

				long startX = getView().bpToTrack(startBp);
				long endX = getView().bpToTrack(endBp);

				drawables.add(new RectDrawable(new Rectangle((int)startX, 0,
						(int)(endX-startX), this.getHeight()), color, color));
			}
		}


		return drawables;
	}

	@Override
	public void processAreaResult(AreaResult areaResult) {

		this.regions.addAll(areaResult.getContents());
		getView().redraw();
	}

	@Override
	public boolean isConcised() {

		return false;
	}

	@Override
	public boolean isStretchable() {

		return false;
	}

	@Override
	public Integer getHeight() {
		if (isVisible()) {
			return 5;
		} else {
			return 0;
		}
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
		datas.put(file, new HashSet<ColumnType>(Arrays.asList(new ColumnType[] {})));
		return datas;
	}
	
    @Override
    public String getName() {
    	return "RepeatMaskerTrack";
    }
}
