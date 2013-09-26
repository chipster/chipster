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

public class PeakSelectable extends Selectable {
	
	private static final int MIN_VISIBLE_PEAK_SIZE = 5;
	public static final int PEAK_SYMBOL_HEIGHT = 16;


	private RectDrawable rectDrawable;	
	private Region region;
	private SelectionText selectionText;

	public PeakSelectable(Region region, IndexKey indexKey, SelectionText selectionText) {		
		super(indexKey);
		this.region = region;
		this.selectionText = selectionText;
	}

	@Override
	public List<Drawable> getDrawables() {
		List<Drawable> drawables = new LinkedList<>();

		if (isSelected()) {

			rectDrawable.color = rectDrawable.color.brighter();
		}
		drawables.add(rectDrawable);
		return drawables;
	}
	
	public Region getRegion() {
		return region;
	}

	public void render(GBrowserView view, Color color) {		

		Rectangle rect = new Rectangle();

		rect.x = view.bpToTrack(region.start);
		rect.width = view.bpToTrack(region.end) - rect.x;

		if (rect.width < MIN_VISIBLE_PEAK_SIZE) {
			rect.width = MIN_VISIBLE_PEAK_SIZE;
		}

		rect.y = 0;
		rect.height = PEAK_SYMBOL_HEIGHT;

		this.rectDrawable = new RectDrawable(rect, color, color.darker());
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