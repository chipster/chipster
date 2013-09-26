package fi.csc.microarray.client.visualisation.methods.gbrowser.track;

import java.awt.Color;
import java.util.LinkedList;

import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserConstants;
import fi.csc.microarray.client.visualisation.methods.gbrowser.gui.GBrowserView;
import fi.csc.microarray.client.visualisation.methods.gbrowser.runtimeIndex.CnaConversion;

public class CnaLogRatiosTrackGroup extends CnaTrackGroup {

	public CnaLogRatiosTrackGroup(GBrowserView view, CnaConversion dataThread,
			LinkedList<String> sampleNames, String title) {
		super(view, dataThread, sampleNames, title + " log ratios", true);
	}
	
	@Override
	public void addTracks() {			
		super.addTracks();
				
		if (!isMinimized()) {
			for (int i = 0; i < sampleNames.size(); i++) {			
				String name = sampleNames.get(i);			
				addLogRatioTrack(name, i);							
			}
		}
	}
	


	private void addLogRatioTrack(String name, int i) {
		TitleTrack title = new TitleTrack(name, Color.black, GBrowserConstants.SCATTERPLOT_TITLE_COLOR);
		title.setView(view);
		addTrack(title);			
		
		ScatterplotTrack logRatio = new ScatterplotTrack(GBrowserConstants.BED_COLOR, 24, -1.0f, 1.0f, i, 0, Long.MAX_VALUE);				
		logRatio.setView(view);
		logRatio.addDataThread(dataThread);		
		addTrack(logRatio);
	}
}
