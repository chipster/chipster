package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.awt.Point;
import java.awt.Rectangle;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.RectDrawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.IndexKey;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.Region;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.SelectionText;

public class ScatterplotPoint extends Selectable {

	private RectDrawable rectDrawable;
	private Region region;
	private Float value;
	private Color color;
	private SelectionText selectionText;

	public ScatterplotPoint(Region region, IndexKey indexKey, Float value, Color color, SelectionText selectionText) {		
		super(indexKey);
		this.region = region;
		this.value = value;
		this.color = color;
		this.selectionText = selectionText;
	}

	@Override
	public List<Drawable> getDrawables() {
		List<Drawable> drawables = new LinkedList<>();
		drawables.add(rectDrawable);
		
		if (isSelected()) {

			Rectangle selectedBounds = rectDrawable.getBounds();
			selectedBounds.grow(4, 4);
			RectDrawable selected = new RectDrawable(selectedBounds, null, Color.black);
			drawables.add(selected);
		}
		return drawables;
	}

	public Region getRegion() {
		return region;
	}

	public void render(GBrowserView view, int minSize, ScatterplotTrack scatterplotTrack) {
			
		int x1 = view.bpToTrack(region.start);
		int width = view.bpToTrack(region.end) - x1;

		width = Math.max(width, minSize);

		int y = scatterplotTrack.getY(value);

		Rectangle rect = new Rectangle();

		rect.y = y;
		rect.x = x1;
		rect.width = width;
		rect.height = minSize;

		this.rectDrawable =  new RectDrawable(rect, color, color);	
	}

	@Override
	public String getText() {
		if (selectionText != null) {
			return selectionText.getText();
		}
		return null;		
	}

	@Override
	public double getDistance(Point point) {
		return this.distance(rectDrawable.getBounds(), point);
	}
}
