package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Point;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.IndexKey;

/**
 * PassiveItem is drawn but doesn't allow user to select it. It can be drawn like any other Selectable object,
 * but it won't get selected when it's clicked. This is useful for track scales, separators and decorations, 
 * because usually there is no need for user to select them. 
 * 
 * @author klemela
 */
public class PassiveItem extends Selectable {

	private List<Drawable> drawables;

	public PassiveItem(List<Drawable> drawables) {
		super(null);
		this.drawables = drawables;
	}

	@Override
	public List<Drawable> getDrawables() {
		return drawables;
	}

	@Override
	public void setSelected(boolean selected) {
	}

	@Override
	public boolean isSelected() {
		return false;
	}

	@Override
	public IndexKey getIndexKey() {
		return null;
	}

	@Override
	public String getText() {
		return null;
	}

	@Override
	public double getDistance(Point point) {
		return Float.MAX_VALUE;
	}
}
