package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Point;
import java.awt.Rectangle;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.Drawable;
import fi.csc.microarray.client.visualisation.methods.gbrowser.message.IndexKey;

public abstract class Selectable {
	
	private IndexKey indexKey;
	private boolean isSelected;

	public Selectable(IndexKey indexKey) {
		this.indexKey = indexKey;
	}

	public abstract List<Drawable> getDrawables();

	public void setSelected(boolean selected) {
		this.isSelected = selected;
	}
	
	public boolean isSelected() {
		return isSelected;
	}

	public IndexKey getIndexKey() {
		return indexKey;
	}

	public abstract String getText();

	public abstract double getDistance(Point point);
	
	public double distance(Rectangle rect, Point p) {
		
		double dx1 = Math.max(rect.getMinX() - p.x, 0);
		double dy1 = Math.max(rect.getMinY() - p.y, 0);
		double dx2 = Math.max(p.x - rect.getMaxX(), 0);
		double dy2 = Math.max(p.y - rect.getMaxY(), 0);
		
		double dx = Math.max(dx1, dx2);
		double dy = Math.max(dy1, dy2);
		
		 return Math.sqrt(dx*dx + dy*dy);
	}

	@Override
	public int hashCode() {
		final int prime = 31;
		int result = 1;
		result = prime * result
				+ ((indexKey == null) ? 0 : indexKey.hashCode());
		return result;
	}

	@Override
	public boolean equals(Object obj) {
		if (this == obj) {
			return true;
		}
		if (obj == null) {
			return false;
		}
		if (!(obj instanceof Selectable)) {
			return false;
		}
		Selectable other = (Selectable) obj;
		if (indexKey == null) {
			if (other.indexKey != null) {
				return false;
			}
		} else if (!indexKey.equals(other.indexKey)) {
			return false;
		}
		return true;
	}
}
