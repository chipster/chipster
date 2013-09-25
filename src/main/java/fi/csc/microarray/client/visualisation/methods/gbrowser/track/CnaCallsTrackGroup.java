package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.CnaConversion;

public class CnaCallsTrackGroup extends CnaTrackGroup {

	public CnaCallsTrackGroup(GBrowserView view, CnaConversion dataThread,
			LinkedList<String> sampleNames, String title) {
		super(view, dataThread, sampleNames, title + " calls", false);
	}
	
	@Override
	public void addTracks() {			
		super.addTracks();
				
		if (!isMinimized()) {
			for (int i = 0; i < sampleNames.size(); i++) {

				String name = sampleNames.get(i);			
				addTracksForName(name, i);										
			}
		}
	}
	
	private void addTracksForName(String name, int i) {
		TitleTrack title = new TitleTrack(name, Color.black);
		title.setView(view);
		addTrack(title);
		
		CnaFlagTrack flag = new CnaFlagTrack(GBrowserConstants.BED_COLOR, i, Color.RED);
		flag.setView(view);
		flag.addDataThread(dataThread);
		addTrack(flag);

		SeparatorTrack separator1 = new SeparatorTrack(Color.gray, 1);
		separator1.setView(view);
		addTrack(separator1);
	}
}
