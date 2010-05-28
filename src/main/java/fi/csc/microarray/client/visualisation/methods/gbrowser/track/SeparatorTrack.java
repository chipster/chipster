package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Arrays;
import java.util.Collection;

import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.LineDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;

public class SeparatorTrack extends Track {

	private Color color;
	private int thickness;
	private long maxBpLength;
	private long minBpLength;

	public SeparatorTrack(View view, Color color, int thickness, long minBpLength, long maxBpLength) {
		super(view, null);
		this.color = color;
		this.thickness = thickness;
		this.minBpLength = minBpLength;
		this.maxBpLength = maxBpLength;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();
		for (int i = 0; i < thickness; i++) {
			drawables.add(new LineDrawable(0, 1+i, getView().getWidth(), 1+i, color));
		}

		return drawables;
	}

	public void processAreaResult(AreaResult areaResult) {
		// ignore
	}

	@Override
	public int getMaxHeight() {

		// the track is hidden if outside given bp length boundaries
		if (getView().getBpRegion().getLength() > minBpLength
				&& getView().getBpRegion().getLength() <= maxBpLength) {

			return 3;

		} else {
			return 0;
		}
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