package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.Arrays;
import java.util.Collection;

import fi.csc.microarray.client.visualisation.methods.gbrowser.View;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.drawable.TextDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.fileFormat.ColumnType;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.AreaResult;

public class TitleTrack extends Track {

	private Color color;
	private String title;

	public TitleTrack(View view, String title, Color color) {
		super(view, null);
		this.color = color;
		this.title = title;
	}

	@Override
	public Collection<Drawable> getDrawables() {
		Collection<Drawable> drawables = getEmptyDrawCollection();
		drawables.add(new TextDrawable(5, 13, title, color));
		return drawables;
	}

	public void processAreaResult(AreaResult areaResult) {
		// ignore
	}

	@Override
	public int getMaxHeight() {
		return 10;
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