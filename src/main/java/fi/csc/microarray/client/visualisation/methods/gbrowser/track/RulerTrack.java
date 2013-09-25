package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Collection;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoord;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.BpCoordDouble;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.DataResult;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;

/**
 * Ruler that shows coordinates.
 *
 */
public class RulerTrack extends Track {

	private static final int textY = 10;
	private final static int boxY = 12;
	private final static int MINOR_STEPS = 10;
	private static final int TICK_W = 2; 
	private static final int BOX_HEIGHT = 4;
	private static final int TICK_H = 4; //on both sides of the box

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

			drawables.add(new TextDrawable(x + 4, textY, text, Color.black));
			//drawables.add(new LineDrawable(x, textY - 2, x, textY + 6, Color.black));
			drawables.add(new RectDrawable(x - TICK_W, boxY - TICK_H, TICK_W, TICK_H * 2 + BOX_HEIGHT, Color.black, Color.black));
		}
		
		boolean whiteStart = false;

		drawables.addAll(getRuler(new Region(start, end, region.start.chr), steps * MINOR_STEPS, whiteStart));

		return drawables;
	}

	private Collection<Drawable> getRuler(Region bpRegion, int steps, boolean whiteStart) {

		Collection<Drawable> drawables = getEmptyDrawCollection();

		boolean isWhite = whiteStart;

		double increment = bpRegion.getLength() / (double) steps;
		BpCoordDouble boxBp = new BpCoordDouble(bpRegion.start);
		
        float lastBoxX = getView().bpToTrackFloat(bpRegion.start);
		float boxX = 0;

		for (int i = 0; i < steps; i++) {

			Color c = (isWhite = !isWhite) ? Color.white : Color.black;
			boxBp = boxBp.move(increment);

			boxX = getView().bpToTrackFloat(boxBp.asBpCoord());

			drawables.add(new RectDrawable(Math.round(lastBoxX), boxY,
			        Math.round(boxX) - Math.round(lastBoxX), BOX_HEIGHT, c, null));

			lastBoxX = boxX;
		}
		
		// Draw a single border
		float startX = getView().bpToTrackFloat(bpRegion.start);
		float endX = boxX;
        drawables.add(new RectDrawable(
                Math.round(startX), boxY,
                Math.round(endX) - Math.round(startX), BOX_HEIGHT, null, Color.black));		
		
		return drawables;

	}

	public void processDataResult(DataResult dataResult) {
		// no data
	}
	
    @Override
    public int getTrackHeight() {
        return boxY + BOX_HEIGHT + TICK_H;
    }
	
	@Override
	public String getTrackName() {
		return "ruler";
	}
}
