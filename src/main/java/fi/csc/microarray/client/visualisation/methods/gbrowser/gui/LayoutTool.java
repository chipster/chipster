package fi.csc.microarray.client.visualisation.methods.gbrowser.gui;

import java.util.Collection;
import java.util.List;

import fi.csc.microarray.client.visualisation.methods.gbrowser.track.Track;
import fi.csc.microarray.client.visualisation.methods.gbrowser.track.TrackGroup;

public class LayoutTool {
	
	public static enum LayoutMode { FIXED, FILL, FULL}

	public static LayoutMode inferTrackGroupLayoutMode(List<Track> tracks) {
		for (Track component : tracks) {
			if (component.isVisible() && component.isSuitableViewLength() && component.getLayoutMode() == LayoutMode.FILL) {
				return LayoutMode.FILL;
			}
		}
		return LayoutMode.FIXED;
	}

	public static LayoutMode inferScrollGroupLayoutMode(List<TrackGroup> trackGroups) {
		for (TrackGroup component : trackGroups) {
			if (component.getLayoutMode() == LayoutMode.FILL) {
				return LayoutMode.FILL;
			}
		}
		return LayoutMode.FIXED;
	}

	public static LayoutMode inferViewLayoutMode(
			Collection<ScrollGroup> scrollGroups) {
		for (ScrollGroup component : scrollGroups) {
			if (component.getLayoutMode() == LayoutMode.FILL) {
				return LayoutMode.FILL;
			}
		}
		return LayoutMode.FIXED;
	}
}
