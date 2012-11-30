package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.util.Collection;
import java.util.LinkedList;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;

public class ScrollGroup implements LayoutComponent, LayoutContainer {

	public List<TrackGroup> trackGroups = new LinkedList<TrackGroup>();
	private int layoutHeight;
	private ScrollPosition defaultScrollPosition = ScrollPosition.MID;
	private String name;
	private int maxViewPortHeight = 300;
	private boolean scrollEnabled = false;;
	
	public ScrollGroup(String name) {
		this.name = name;
	}
	
	public ScrollGroup(String name, ScrollPosition defaultScrollPosition) {
		this(name);
		this.defaultScrollPosition = defaultScrollPosition;
	}
	
	public ScrollGroup(String name, ScrollPosition defaultScrollPosition, int maxViewPortHeight) {
		this(name, defaultScrollPosition);
		this.maxViewPortHeight = maxViewPortHeight;
	}

	public ScrollGroup() {
	}

	public ScrollGroup(String name, boolean scrollEnabled) {
		this(name);
		this.scrollEnabled = scrollEnabled;
	}

	@Override
	public Collection<? extends LayoutComponent> getLayoutComponents() {
		return trackGroups;
	}

	@Override
	public boolean isFixedHeight() {
		return LayoutTool.isFixedHeight(this);
	}

	@Override
	public int getHeight() {
		return layoutHeight;
	}

	@Override
	public void setHeight(int height) {
		this.layoutHeight = height;
	}

	@Override
	public int getMinHeight() {
		return 0;
	}

	@Override
	public boolean isVisible() {
		return true;
	}

	public void addTrack(Track track) {
		trackGroups.add(new TrackGroup(track));
	}

	public void addTrackGroup(TrackGroup group) {
		trackGroups.add(group);
	}

	public Collection<TrackGroup> getTrackGroups() {
		return trackGroups;
	}

	public int getFullHeight() {
		return 0;
	}

	public ScrollPosition getDefaultScrollPosition() {
		return defaultScrollPosition ;
	}
	
	public enum ScrollPosition { START, MID };
	
	public String toString() {
		return ScrollGroup.class + " " + name;
	}

	//FIXME can be removed?
	public int getMaxViewPortHeight() {
		return maxViewPortHeight ;
	}

	public boolean isScrollEnabled() {
		return scrollEnabled;
	}

	@Override
	public int getCanvasHeight() {
		return LayoutTool.getCanvasHeight(this);
	}

	public int getScrollReferenceY() {
		//FIXME create generic solution, for example extend this class for different special cases (samples, annotations, etc. )
		if ("Annotations".equals(name)) {
			return Math.max(trackGroups.get(0).getTracks().get(1).getCanvasHeight(), trackGroups.get(0).getTracks().get(1).getHeight());
		}
		return 0;
	}
}
