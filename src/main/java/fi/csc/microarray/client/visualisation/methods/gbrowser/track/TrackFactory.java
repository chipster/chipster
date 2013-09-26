package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserPlot;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;

/**
 * Utility class for creating predefined {@link TrackGroup} objects.  
 *
 *@author Petri Klemelä, Aleksi Kallio
 *
 */
public class TrackFactory {	
    
	public static TrackGroup getThinSeparatorTrackGroup(GBrowserPlot genomePlot) {
		GBrowserView view = genomePlot.getDataView();
		TrackGroup group = new TrackGroup(view);
		SeparatorTrack separator = new SeparatorTrack(Color.LIGHT_GRAY, 3); 
		separator.setView(view);
		group.addTrack(separator);
		return group;
	}

	public static TrackGroup getThickSeparatorTrackGroup(GBrowserPlot genomePlot) {
		GBrowserView view = genomePlot.getDataView();
		TrackGroup group = new TrackGroup(view);
		
		SeparatorTrack3D separator1 = new SeparatorTrack3D(false);		
		separator1.setView(view);	
		group.addTrack(separator1);
		
		EmptyTrack empty = new EmptyTrack(2);
		empty.setView(view);		
		group.addTrack(empty);
		
		SeparatorTrack3D separator2 = new SeparatorTrack3D(true);		
		separator2.setView(view);	
		group.addTrack(separator2);
		
		return group;
	}
}
