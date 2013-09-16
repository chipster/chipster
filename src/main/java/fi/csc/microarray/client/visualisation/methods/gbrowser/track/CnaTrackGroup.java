package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.CnaConversion;

public abstract class CnaTrackGroup extends TrackGroup {

	protected CnaConversion dataThread;
	protected LinkedList<String> sampleNames;

	public CnaTrackGroup(GBrowserView view, CnaConversion dataThread,
			LinkedList<String> sampleNames, String title, boolean minimized) {
		super(view);
		
		this.dataThread = dataThread;
		this.sampleNames = sampleNames;
		
		setMinimzed(minimized);
		this.setName(title);
		
		getStatusAnimation().addDataThread(dataThread);
		setSettingsEnabled(true);		
		addTracks();
	}
	
	@Override
	public void addTracks() {
		tracks.clear();
		SeparatorTrack separator = new SeparatorTrack(Color.white, 20);
		separator.setView(view);
		addTrack(separator);
	}
}
