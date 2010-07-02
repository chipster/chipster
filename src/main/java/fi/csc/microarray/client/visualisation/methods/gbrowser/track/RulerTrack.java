package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordRegion;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.RegionContent;

/**
 * Ruler that shows coordinates.
 *
 */
public class RulerTrack extends Track {

	private static final int textY = 10;
	private final static int MINOR_STEPS = 10;

	List<Long> info = new ArrayList<Long>();

	public RulerTrack(View view) {
		super(view, null);
	}

	@Override
	public Collection<Drawable> getDrawables() {

		Collection<Drawable> drawables = getEmptyDrawCollection();
		BpCoordRegion region = getView().getBpRegion();

		long magnitude = (long) Math.pow(10, (int) Math.log10(region.getLength()));

		final long start = region.start.bp - region.start.bp % magnitude;
		final int steps = (int) Math.ceil(region.getLength() / magnitude) + 1;
		final long end = start + steps * magnitude;

		for (long bp = start; bp <= end; bp += magnitude) {

			int x = getView().bpToTrack(new BpCoord(bp, region.start.chr));
			String text = Utils.toHumanReadable(bp);

			drawables.add(new TextDrawable(x, textY, text, Color.black));
		}

		drawables.addAll(getRuler(new BpCoordRegion(start, end, region.start.chr), steps * MINOR_STEPS, start / magnitude % 2 == 1));

		return drawables;
	}

	private Collection<Drawable> getRuler(BpCoordRegion bpRegion, int steps, boolean whiteStart) {

		Collection<Drawable> drawables = getEmptyDrawCollection();

		boolean isWhite = whiteStart;
		final int boxHeight = 5;

		double increment = bpRegion.getLength() / (double) steps;
		BpCoordDouble boxBp = new BpCoordDouble(bpRegion.start);
		int boxX;
		int lastBoxX = getView().bpToTrack(bpRegion.start);

		info.clear();

		for (int i = 0; i < steps; i++) {

			Color c = (isWhite = !isWhite) ? Color.white : Color.black;
			boxBp = boxBp.move(increment);

			info.add((long) (double) boxBp.bp);

			boxX = getView().bpToTrack(boxBp.asBpCoord());

			drawables.add(new RectDrawable(lastBoxX, textY, boxX - lastBoxX, boxHeight, c, Color.black));

			Color lineColor = (i % MINOR_STEPS == MINOR_STEPS - 1) ? new Color(0, 0, 0, 64) : new Color(0, 0, 0, 32);

			drawables.add(new LineDrawable(boxX, -getView().getHeight() + getHeight(), boxX, getHeight(), lineColor));

			lastBoxX = boxX;
		}
		return drawables;

	}

	public void processAreaResult(AreaResult<RegionContent> areaResult) {
		// no data
	}

	public List<Long> getRulerInfo() {
		return info;
	}
	
    @Override
    public Integer getHeight() {
        return textY * 2;
    }

	@Override
	public boolean isStretchable () {
		return false;
	}

	@Override
	public Collection<ColumnType> getDefaultContents() {
		return Arrays.asList(new ColumnType[] {});
	}

	@Override
	public boolean isConcised() {
		return false;
	}
}
