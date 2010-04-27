package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.io.File;
import java.util.Arrays;
import java.util.Collection;
import java.util.Iterator;
import java.util.TreeSet;

import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.FileParser;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

public class ProfileTrack extends Track {

	private Collection<RegionContent> peaks = new TreeSet<RegionContent>();

	private long maxBpLength;
	private long minBpLength;

	private boolean wasLastConcised = true;
	private Color color;


	public ProfileTrack(View view, File file, Class<? extends AreaRequestHandler> handler, FileParser inputParser, Color color, long minBpLength, long maxBpLength) {
		super(view, file, handler, inputParser);
		this.color = color;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();

		if (peaks != null) {

			Iterator<RegionContent> iter = peaks.iterator();
			RegionContent last = null;
			while (iter.hasNext()) {

				RegionContent value = iter.next();

				if (!value.region.intercepts(getView().getBpRegion())) {
					iter.remove();
					continue;
				}
				
				if (last != null) {
					createDrawable(value.region.start, value.region.start.bp.intValue() % 10, last.region.start, last.region.start.bp.intValue() % 10, color, drawables);
				}
				last = value;
			}
		}

		return drawables;
	}

	private void createDrawable(BpCoord startBp, int h, BpCoord endBp, int h2, Color c, Collection<Drawable> drawables) {

		int x = getView().bpToTrack(startBp);
		int x2 = getView().bpToTrack(endBp);

		int y = (int) (getView().getTrackHeight()) + h;
		int y2 = (int) (getView().getTrackHeight()) + h2;

		drawables.add(new LineDrawable(x, y, x2, y2, c));
	}

	public void processAreaResult(AreaResult<RegionContent> areaResult) {

		if (areaResult.status.concise == this.isConcised()) {
			this.peaks.add(areaResult.content);
			getView().redraw();
		}
	}


	@Override
	public void updateData() {
		if (wasLastConcised != isConcised()) {
			peaks.clear();
			wasLastConcised = isConcised();
		}
		super.updateData();
	}

	@Override
	public int getMaxHeight() {
		if (getView().getBpRegion().getLength() > minBpLength && getView().getBpRegion().getLength() <= maxBpLength) {
			return super.getMaxHeight();
			
		} else {
			return 0;
		}
	}

	@Override
	public Collection<ColumnType> getDefaultContents() {
		return Arrays.asList(new ColumnType[] { });
	}

	@Override
	public boolean isConcised() {
		return false;
	}
}
