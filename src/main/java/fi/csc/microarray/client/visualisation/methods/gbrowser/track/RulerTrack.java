package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;
import java.util.Map;
import java.util.Set;

import fi.csc.microarray.client.visualisation.methods.gbrowser.dataFetcher.AreaRequestHandler;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

/**
 * Ruler that shows coordinates.
 *
 */
public class RulerTrack extends Track {

	private static final int textY = 10;
	private final static int MINOR_STEPS = 10;

	public RulerTrack() {
	}

	@Override
	public Collection<Drawable> getDrawables() {

		Collection<Drawable> drawables = getEmptyDrawCollection();
		Region region = getView().getBpRegion();

		long magnitude = (long) Math.pow(10, (int) Math.log10(region.getLength()));

		final long start = region.start.bp - region.start.bp % magnitude;
		final int steps = (int) Math.ceil(region.getLength() / (double)magnitude) + 1;
		final long end = start + steps * magnitude;

		for (long bp = start; bp <= end; bp += magnitude) {

			int x = getView().bpToTrack(new BpCoord(bp, region.start.chr));
			String text = Utils.toHumanReadable(bp);

			drawables.add(new TextDrawable(x, textY, text, Color.black));
		}
		
		boolean whiteStart = false;

		drawables.addAll(getRuler(new Region(start, end, region.start.chr), steps * MINOR_STEPS, whiteStart));

		return drawables;
	}

	private Collection<Drawable> getRuler(Region bpRegion, int steps, boolean whiteStart) {

		Collection<Drawable> drawables = getEmptyDrawCollection();

		boolean isWhite = whiteStart;
		final int boxHeight = 5;

		double increment = bpRegion.getLength() / (double) steps;
		BpCoordDouble boxBp = new BpCoordDouble(bpRegion.start);
		
        float lastBoxX = getView().bpToTrackFloat(bpRegion.start);
		float boxX = 0;

		for (int i = 0; i < steps; i++) {

			Color c = (isWhite = !isWhite) ? Color.white : Color.black;
			boxBp = boxBp.move(increment);

			boxX = getView().bpToTrackFloat(boxBp.asBpCoord());

			drawables.add(new RectDrawable(Math.round(lastBoxX), textY,
			        Math.round(boxX) - Math.round(lastBoxX), boxHeight, c, null));

			lastBoxX = boxX;
		}
		
		// Draw a single border
		float startX = getView().bpToTrackFloat(bpRegion.start);
		float endX = boxX;
        drawables.add(new RectDrawable(
                Math.round(startX), textY,
                Math.round(endX) - Math.round(startX), boxHeight, null, Color.black));		
		
		return drawables;

	}

	public void processAreaResult(AreaResult areaResult) {
		// no data
	}
	
    @Override
    public int getHeight() {
        return textY * 2;
    }

    @Override
    public Map<AreaRequestHandler, Set<ColumnType>> requestedData() {
        return null;
    }
	
	@Override
	public String getName() {
		return "ruler";
	}
}
