package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.DataThread;

public class RegionTrackGroup extends TrackGroup {
	
	private RegionTrack regionTrack;
	private DataThread dataThread;

	public RegionTrackGroup(GBrowserView view, DataThread dataThread, String title) {
		super(view);
		this.setSettingsEnabled(true);
		this.dataThread = dataThread;
		
		this.setName(title);
		addTracks();	
	}
	
	@Override
	public void addTracks() {		
		tracks.clear();

		getStatusAnimation().addDataThread(dataThread);

		if (!isMinimized()) {
			SeparatorTrack separator = new SeparatorTrack(Color.white, 20);
			separator.setView(view);
			addTrack(separator);
			regionTrack = new RegionTrack(GBrowserConstants.BED_COLOR);
			regionTrack.setView(view);
			regionTrack.addDataThread(dataThread);
			addTrack(regionTrack);			
		}
	}
}
